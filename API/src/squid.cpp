#include "squid.h"

using namespace std;

namespace constants
{
    // DEBUG PARAMETERS
    const int DEBUG = 0;
    
    // BGV PARAMETERS
    // Plaintext prime modulus
    const unsigned long P = 131;
    // Cyclotomic polynomial - defines phi(m)
    const unsigned long M = 130;//17293;
    // Hensel lifting (default = 1)
    const unsigned long R = 1;
    // Number of bits of the modulus chain
    const unsigned long BITS = 1000;//431;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    const unsigned long C = 3;
}

template<typename T, typename Allocator>
void print_vector(const vector<T, Allocator>& vect, int num_entries)
{
    cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++){
        cout << ", " << vect[i];
    }
    cout << endl;
}

// Function for extended Euclidean Algorithm
int gcdExtended(int a, int b, int* x, int* y){
 
    // Base Case
    if (a == 0) {
        *x = 0, *y = 1;
        return b;
    }
 
    // To store results of recursive call
    int x1, y1;
    int gcd = gcdExtended(b % a, a, &x1, &y1);
 
    // Update x and y using results of recursive
    // call
    *x = y1 - (b / a) * x1;
    *y = x1;
 
    return gcd;
}
// Function to find modulo inverse of a
int modInverse(int A, int M){
    int x, y;
    int g = gcdExtended(A, M, &x, &y);
    if (g != 1)
        throw invalid_argument("No inverse");
    else {
 
        // m is added to handle negative x
        int res = (x % M + M) % M;
        return res;
    }
}


int get_inverse(int nom, int dom, int p){
    int bottom_inverse = (modInverse(dom, (int)p)) % (int)p;
    if (bottom_inverse < 0){
        bottom_inverse += (int)p;
    }
    int result = (nom * bottom_inverse) % p;
    if (result < 0){
        result += p;
    }
    return result;
}

void Squid::SetServerToExample(){
    vector<vector<unsigned long>> fake_db = vector<vector<unsigned long>> {
          {0, 0, 0, 1, 1, 1, 2, 2, 2, 0},
          {0, 1, 2, 0, 1, 2, 0, 1, 2, 1},
          {0, 0, 0, 0, 1, 1, 0, 1, 0, 0}
     };
    vector<string> headers = vector<string>{"snp1", "snp2","ALS"};
    SetColumnHeaders(headers);
    cout << "set headers" << endl;
    SetData(fake_db);
}

Squid::Squid(): context(helib::ContextBuilder<helib::BGV>()
                               .m(constants::M)
                               .p(constants::P)
                               .r(constants::R)
                               .bits(constants::BITS)
                               .c(constants::C)
                               .build()), secret_key(context){
    cout << "called" << endl;

    secret_key.GenSecKey();

    helib::addSome1DMatrices(secret_key);

    const helib::EncryptedArray& ea = context.getEA();
    num_slots = ea.size();
    plaintext_modulus = context.getP();

    one_over_two = get_inverse(1,2,plaintext_modulus);
    neg_three_over_two = get_inverse(-3,2,plaintext_modulus);
    neg_one_over_two = get_inverse(-1,2,plaintext_modulus);

    db_set = false;
    helib::Ptxt<helib::BGV> ptxt(context);
    public_key_ptr = new helib::PubKey(secret_key);

    SetServerToExample();
    PrintContext();
}

void Squid::GenData(int _num_rows, int _num_cols){
    num_rows = _num_rows;
    num_cols = _num_cols;

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    encrypted_db = vector<vector<helib::Ctxt>>();
    for(int i = 0; i < num_cols; i++){
        vector<helib::Ctxt> cipher_vector = vector<helib::Ctxt>();
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ptxt<helib::BGV> ptxt(context);

            helib::Ctxt ctxt(*public_key_ptr);

            public_key_ptr->Encrypt(ctxt, ptxt);

            cipher_vector.push_back(ctxt);
        }
        encrypted_db.push_back(cipher_vector);
    }

    db_set = true;
}

const helib::Context& Squid::GetContext() const{
    return context;
}

int Squid::GetNumRows() const{
    return num_rows;
}

const helib::SecKey& Squid::GetSK() const{
    return secret_key;
}




const vector<string>& Squid::GetColumnHeaders() const{
    return column_headers;
}


void Squid::SetData(vector<vector<unsigned long>> &db){
    num_cols = db.size();
    if (num_cols == 0){
        throw invalid_argument("ERROR: DB has zero columns! THIS DOES NOT WORK!");
    }
    cout << __LINE__ << endl;
    num_rows = db[0].size();

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    encrypted_db = vector<vector<helib::Ctxt>>();
    
    for(int i = 0; i < num_cols; i++){
            cout << __LINE__ << endl;

        vector<helib::Ctxt> cipher_vector = vector<helib::Ctxt>();
        for (int j = 0; j < num_compressed_rows; j++){
                cout << __LINE__ << endl;


            helib::Ptxt<helib::BGV> ptxt(context);
                cout << __LINE__ << endl;


            int entries_left = min(num_slots, num_rows - (j * num_slots));
            for (int k = 0; k < entries_left; k++){
                ptxt[k] = db[i][j*num_slots + k];
            }

            helib::Ctxt ctxt(*public_key_ptr);


            public_key_ptr->Encrypt(ctxt, ptxt);

            cipher_vector.push_back(ctxt);
        }
        encrypted_db.push_back(cipher_vector);
    }
    db_set = true;
}

void Squid::SetColumnHeaders(vector<string> &headers){
    column_headers = vector<string>();
    for (int i = 0; i < headers.size(); i++){
        column_headers.push_back(headers[i]);
    }
}

helib::Ctxt Squid::CountingQuery(bool conjunctive, vector<pair<int, int>>& query) const{
    if (!db_set){
        throw invalid_argument("ERROR: DB needs to be set to run query");
    }
    vector<vector<helib::Ctxt>> cols = filter(query);

    int num_columns = cols[0].size();

    vector<helib::Ctxt> filter_results;
    if (conjunctive){
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else{
        for (int i = 0; i < num_compressed_rows; i++){
            for (int j = 0; j < num_columns; j++){
                AddOneMod2(cols[i][j]);
            }
        }
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (int j = 0; j < num_compressed_rows; j++){
            AddOneMod2(filter_results[j]);
        }
    }

    helib::Ctxt result = AddMany(filter_results);
    print_vector(Decrypt(result), 10);

    //print_vector(Decrypt(result), num_slots);
    result = SquashCtxtWithMask(result, 0);

    return result;


}

pair<helib::Ctxt, helib::Ctxt> Squid::MAFQuery(int snp, bool conjunctive, vector<pair<int, int>> &query) const{
    vector<vector<helib::Ctxt>> cols = filter(query);
    int num_columns = cols[0].size();

    vector<helib::Ctxt> filter_results;
    if (conjunctive){
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else{
        for (int i = 0; i < num_compressed_rows; i++){
            for (int j = 0; j < num_columns; j++){
                AddOneMod2(cols[i][j]);
            }
        }
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (int j = 0; j < num_compressed_rows; j++){
            AddOneMod2(filter_results[j]);
        }
    }

    vector<helib::Ctxt> indv_MAF = vector<helib::Ctxt>();

    for (int i = 0; i < num_compressed_rows; i++){
        helib::Ctxt clone = encrypted_db[snp][i];
        clone *= filter_results[i];
        indv_MAF.push_back(clone);
    }

    helib::Ctxt freq = AddMany(indv_MAF);
    helib::Ctxt number_of_patients = AddMany(filter_results);

    freq = SquashCtxtWithMask(freq, 0);
    number_of_patients = SquashCtxtWithMask(number_of_patients, 0);

    number_of_patients.multByConstant(NTL::ZZX(2));
    return pair(freq, number_of_patients);
}

vector<helib::Ctxt> Squid::PRSQuery(vector<pair<int, int>>& prs_params) const{

    vector<helib::Ctxt> scores;

    for(int j = 0; j < num_compressed_rows; j++){
        vector<helib::Ctxt> indvs_scores;
        for(pair<int, int> i : prs_params){
            helib::Ctxt temp = encrypted_db[i.first][j];

            temp.multByConstant(NTL::ZZX(i.second));
            indvs_scores.push_back(temp);
        }
        helib::Ctxt score = AddMany(indvs_scores);
        scores.push_back(score);
    }
    return scores;
}

vector<pair<helib::Ctxt, helib::Ctxt>> Squid::ChiSquareQuery(int disease_column, int number_of_chi){
    vector<pair<helib::Ctxt, helib::Ctxt>> chi_square_results = vector<pair<helib::Ctxt, helib::Ctxt>>();

    helib::Ctxt n11(*public_key_ptr);
    helib::Ctxt c1(*public_key_ptr);

    // r1
    helib::Ctxt r1 = AddManySafe(encrypted_db[disease_column]);
    SquashCtxtWithMask(r1,0);
    CtxtExpand(r1);
    r1.multByConstant(NTL::ZZX(2));
    // r1


    // Y
    for (int c = 0; c < number_of_chi; c++){
        vector<helib::Ctxt> ytS = vector<helib::Ctxt>();
        for (int r = 0; r < num_compressed_rows; r++){
            helib::Ctxt s = encrypted_db[c][r];
            s.multiplyBy(encrypted_db[disease_column][r]);
            ytS.push_back(s);
        }
        helib::Ctxt sum = AddManySafe(ytS);
        n11.addCtxt(SquashCtxtWithMask(sum, c));
        sum = AddManySafe(encrypted_db[c]);
        c1.addCtxt(SquashCtxtWithMask(sum,c));
    }

    // D
    helib::Ctxt d = Encrypt(2 * num_rows);

    //chisquare

    helib::Ctxt rc = r1;
    rc.multiplyBy(c1);

    helib::Ctxt d_mins_c1 = d;
    d_mins_c1 -= c1;
    helib::Ctxt d_mins_r1 = d;
    d_mins_r1 -= r1;

    helib::Ctxt n11_times_d = d;
    n11_times_d *= n11;

    n11_times_d -= rc;
    n11_times_d.square();
    n11_times_d *= d;

    helib::Ctxt num = n11_times_d;

    helib::Ctxt den = rc;
    den *= d_mins_c1;
    den *= d_mins_r1;

    chi_square_results.push_back(pair(num, den));

    return chi_square_results;
}

vector<pair<helib::Ctxt, helib::Ctxt>> Squid::ChiSquareQuery(bool conjunctive, vector<pair<int, int>>& query, int disease_column, int number_of_chi){
    if (!db_set){
        throw invalid_argument("ERROR: DB needs to be set to run query");
    }

    vector<vector<helib::Ctxt>> cols = filter(query);

    int num_columns = cols[0].size();

    vector<helib::Ctxt> filter_results;
    if (conjunctive){
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else{
        for (int i = 0; i < num_compressed_rows; i++){
            for (int j = 0; j < num_columns; j++){
                AddOneMod2(cols[i][j]);
            }
        }
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (int j = 0; j < num_compressed_rows; j++){
            AddOneMod2(filter_results[j]);
        }
    }

    vector<helib::Ctxt> y = vector<helib::Ctxt>();
    vector<helib::Ctxt> un_minus_y = vector<helib::Ctxt>();

    for (int j = 0; j < num_compressed_rows; j++){
        helib::Ctxt predicate = filter_results[j];

        y.push_back(predicate);

        helib::Ctxt inverse_predicate = predicate;

        AddOneMod2(inverse_predicate);
        un_minus_y.push_back(inverse_predicate);
    }


    vector<pair<helib::Ctxt, helib::Ctxt>> chi_square_results = vector<pair<helib::Ctxt, helib::Ctxt>>();

    for (int i = 0; i < num_cols - 1; i++){
        for (int j = 0; j < num_compressed_rows; j++){

        }
    }

    return chi_square_results;
}


void Squid::AddOneMod2(helib::Ctxt& a) const{
    //   0 -> 1
    //   1 -> 0
    // f(x)-> -x+1

    a.negate();
    a.addConstant(NTL::ZZX(1));
}

helib::Ctxt Squid::MultiplyMany(vector<helib::Ctxt>& v) const{
    int num_entries = v.size();
    int depth = ceil(log2(num_entries));

    for (int d = 0; d < depth; d++){
        int jump_factor = pow(2, d);
        int skip_factor = 2 * jump_factor;

        for (int i = 0; i < num_entries; i+= skip_factor){
            v[i].multiplyBy(v[i + jump_factor]);
        }
     }
     return v[0];
}



helib::Ctxt Squid::AddMany(vector<helib::Ctxt>& v) const{
    int num_entries = v.size();
    int depth = ceil(log2(num_entries));

    for (int d = 0; d < depth; d++){
        int jump_factor = pow(2, d);
        int skip_factor = 2 * jump_factor;

        for (int i = 0; i < num_entries; i+= skip_factor){
            v[i] += v[i + jump_factor];
        }
     }
     return v[0];
}

helib::Ctxt Squid::AddManySafe(vector<helib::Ctxt>& v) const{

    helib::Ctxt result(*public_key_ptr);

    for (int i = 0; i < v.size(); i++){
        result += v[i];
    }
    return result;
}

helib::Ctxt Squid::SquashCtxt(helib::Ctxt& ciphertext, int num_data_elements) const{
    const helib::EncryptedArray& ea = context.getEA();

    helib::Ctxt result = ciphertext;

    for (int i = 1; i < num_data_elements; i++) {
        ea.rotate(ciphertext, -(1));
        result += ciphertext;
    }
    return result;
}

helib::Ctxt Squid::SquashCtxtLogTime(helib::Ctxt& ciphertext) const{
    const helib::EncryptedArray& ea = context.getEA();

    int depth = floor(log2(num_slots));

    //cout << "Depth: " <<  depth << endl;

    helib::Ptxt<helib::BGV> mask(context);
    helib::Ptxt<helib::BGV> inverse_mask(context);

    int largest_power_of_two_less_than_or_equal_two_slotsize = 1 << depth;
    for (int i = 0; i < num_slots; i++){
        if (i < largest_power_of_two_less_than_or_equal_two_slotsize){
            mask[i] = 1;
            inverse_mask[i] = 0;
        }
        else{
            mask[i] = 0;
            inverse_mask[i] = 1;
        }
    }
    //cout << "far_end:" << endl;
    helib::Ctxt far_end = ciphertext;
    far_end.multByConstant(inverse_mask);

    //print_vector(Decrypt(far_end), num_slots);

    ea.rotate(far_end, -largest_power_of_two_less_than_or_equal_two_slotsize);

    ciphertext.multByConstant(mask);
    ciphertext += far_end;

    //print_vector(Decrypt(ciphertext), 10);


    for (int d = depth - 1; d >= 0; d--){
        int shift = 1 << d;
        helib::Ctxt clone = ciphertext;
        ea.rotate(clone, (-shift));
        ciphertext += clone;
        //cout << "Printing mid way through with depth:" << d << endl;
        //print_vector(Decrypt(ciphertext), 10);
    }

    return ciphertext;
}

helib::Ctxt Squid::SquashCtxtWithMask(helib::Ctxt& ciphertext, int index) const{
    ciphertext = SquashCtxtLogTime(ciphertext);
    //cout << "Running squash with mask" << endl;
    //print_vector(Decrypt(ciphertext), 10);
    const helib::EncryptedArray& ea = context.getEA();

    ea.rotate(ciphertext, index);

    helib::Ptxt<helib::BGV> mask(context);
    mask[index] = 1;

    ciphertext.multByConstant(mask);

    return ciphertext;
}

void Squid::CtxtExpand(helib::Ctxt &ciphertext) const{
    const helib::EncryptedArray& ea = context.getEA();

    int depth = floor(log2(num_slots));

    helib::Ptxt<helib::BGV> mask(context);

    int largest_power_of_two_less_than_or_equal_two_slotsize = 1 << depth;
    for (int i = 0; i < num_slots - largest_power_of_two_less_than_or_equal_two_slotsize; i++){
        mask[i] = 1;
    }
    for (int d = 0; d < depth; d++){
        int shift = 1 << d;
        helib::Ctxt clone = ciphertext;
        ea.rotate(clone, (shift));
        ciphertext += clone;
    }

    helib::Ctxt clone = ciphertext;
    clone.multByConstant(mask);
    ea.rotate(clone, largest_power_of_two_less_than_or_equal_two_slotsize);
    ciphertext += clone;
}

helib::Ctxt Squid::EQTest(unsigned long a, const helib::Ctxt& b) const{

    helib::Ctxt clone = b;
    helib::Ctxt result = b;

    switch (a){
        case 0:
        {
            //f(x) = x^2 / 2 - 3/2 x + 1
            // 0 -> 1
            // 1 -> 0
            // 2 -> 0
            result.square();

            result.multByConstant(NTL::ZZX(one_over_two));
            clone.multByConstant(NTL::ZZX(neg_three_over_two));

            result += clone;

            result.addConstant(NTL::ZZX(1));
            return result;
        }
        case 1:
        {
            //f(x) = -x^2 + 2x
            // 0 -> 0
            // 1 -> 1
            // 2 -> 0

            clone.square();
            result.multByConstant(NTL::ZZX(2));

            result -= clone;
            return result;
        }
        case 2:{
            //f(x) = x^2 / 2 - x / 2
            // 0 -> 0
            // 1 -> 0
            // 2 -> 1

            result.square();

            result.multByConstant(NTL::ZZX(one_over_two));
            clone.multByConstant(NTL::ZZX(neg_one_over_two));

            result += clone;
            return result;
        }
        default:
            cout << "Can't use a value of a other than 0, 1, or 2" << endl;
            throw invalid_argument("ERROR: invalid value for EQTest");
    }
}

vector<vector<helib::Ctxt>> Squid::filter(vector<pair<int, int>>& query) const{
    vector<vector<helib::Ctxt>> feature_cols;

    for(int j = 0; j < num_compressed_rows; j++){
        vector<helib::Ctxt> indv_vector;
        for(pair<int, int> i : query){

            indv_vector.push_back(EQTest(i.second, encrypted_db[i.first][j]));
        }
        feature_cols.push_back(indv_vector);
    }
    return feature_cols;
}


vector<long> Squid::Decrypt(helib::Ctxt ctxt) const{
    helib::Ptxt<helib::BGV> new_plaintext_result(context);
    secret_key.Decrypt(new_plaintext_result, ctxt);

    vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();

    vector<long> result = vector<long>(num_slots);

    for (int i = 0; i < num_slots; i++){
        result[i] = (long)poly_mod_result[i];
    }

    return result;
}

helib::Ptxt<helib::BGV> Squid::DecryptPlaintext(helib::Ctxt ctxt) const{

    if (constants::DEBUG && ctxt.capacity() < 2){
        cout << "NOISE BOUNDS EXCEEDED!!!" << endl;
    }

    helib::Ptxt<helib::BGV> new_plaintext_result(context);
    secret_key.Decrypt(new_plaintext_result, ctxt);

    return new_plaintext_result;
}

helib::Ctxt Squid::Encrypt(unsigned long a) const{
    helib::Ptxt<helib::BGV> ptxt(context);

    for (int i = 0; i < num_slots; i++)
        ptxt[i] = a;

    helib::Ctxt ctxt(*public_key_ptr);
    public_key_ptr->Encrypt(ctxt, ptxt);

    return ctxt;
}

helib::Ctxt Squid::Encrypt(vector<unsigned long> a) const{
    if (a.size() > num_slots){
        throw invalid_argument("Trying to encrypt vector with too many elements");
    }
    helib::Ptxt<helib::BGV> ptxt(context);

    for (size_t i = 0; i < a.size(); ++i) {
        ptxt[i] = a[i];
    }

    helib::Ctxt ctxt(*public_key_ptr);
    public_key_ptr->Encrypt(ctxt, ptxt);

    return ctxt;
}

helib::Ctxt Squid::GetAnyElement() const{
    return encrypted_db[0][0];
}

void Squid::PrintContext() const{
    context.printout();
    cout << endl;
    cout << "Security: " << context.securityLevel() << endl;
    cout << "Num slots: " << num_slots << endl;
}

string Squid::PrintEncryptedDB(bool with_headers) const{
    string s = "";
        if (with_headers){
        vector<int> string_length_count = vector<int>();

        s = "|";

        for (int i = 0; i < num_cols; i++){
            s += column_headers[i] + "|";
            string_length_count.push_back(column_headers[i].length());
        }
        s += "\n";
        s += "--------------";
        for (int j = 0; j < num_compressed_rows; j++){

            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (int i = 0; i < num_cols; i++){
                temp_storage.push_back(Decrypt(encrypted_db[i][j]));
            }
            for (int jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++){

                s += "\n|";

                for (int i = 0; i < num_cols; i++){
                    if(i>0){
                        for (int space = 0; space < string_length_count[i]; space++){
                            s += " ";
                        }

                    }

                    s += to_string(temp_storage[i][jj]);

                    if(i==num_cols-1){
                        for (int space = 0; space < string_length_count[i]; space++){
                            s += " ";
                        }
                        s += "|";
                    }
                }
            }
        }
    }
    else{
        for (int j = 0; j < num_compressed_rows; j++){

            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (int i = 0; i < num_cols; i++){
                temp_storage.push_back(Decrypt(encrypted_db[i][j]));
            }
            for (int jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++){

                s += "\n|";

                for (int i = 0; i < num_cols; i++){
                    if(i>0){
                        s += " ";
                    }

                    s += to_string(temp_storage[i][jj]);

                    if(i==num_cols-1){
                        s += "|";
                    }
                }
            }
        }
    }
    s += "\n";

    return s;
}

int Squid::GetSlotSize() const{
    return num_slots;
}

// IMPORTED FROM HELIB SOURCE CODE
inline long estimateCtxtSize(const helib::Context& context, long offset)
{
  // Return in bytes.

  // We assume that the size of each element in the DCRT is BINIO_64BIT

  // sizeof(BINIO_EYE_CTXT_BEGIN) = 4;
  // BINIO_32BIT = 4
  // sizeof(long) = BINIO_64BIT = 8
  // xdouble = s * sizeof(long) = 2 * BINIO_64BIT = 16

  // We assume that primeSet after encryption is context.ctxtPrimes
  // We assume we have exactly 2 parts after encryption
  // We assume that the DCRT prime set is the same as the ctxt one

  long size = 0;

  // Header metadata
  size += 24;

  // Begin eye-catcher
  size += 4;

  // Begin Ctxt metadata
  // 64 = header_size = ptxtSpace (long) + intFactor (long) + ptxtMag (xdouble)
  //                    + ratFactor (xdouble) + noiseBound (xdouble)
  size += 64;

  // primeSet.write(str);
  // size of set (long) + each prime (long)
  size += 8 + context.getCtxtPrimes().card() * 8;

  // Begin Ctxt content size
  // write_raw_vector(str, parts);
  // Size of the parts vector (long)
  size += 8;

  long part_size = 0;
  // Begin CtxtPart size

  // skHandle.write(str);
  // powerOfS (long) + powerOfX (long) + secretKeyID (long)
  part_size += 24;
    // Begin DCRT size computation

  // this->DoubleCRT::write(str);
  // map.getIndexSet().write(str);
  // size of set (long) + each prime (long)
  part_size += 8 + context.getCtxtPrimes().card() * 8;

  // DCRT data write as write_ntl_vec_long(str, map[i]);
  // For each prime in the ctxt modulus chain
  //    size of DCRT column (long) + size of each element (long) +
  //    size of all the slots (column in DCRT) (PhiM long elements)
  long dcrt_size = (8 + 8 * context.getPhiM()) * context.getCtxtPrimes().card();

  part_size += dcrt_size;

  // End DCRT size
  // End CtxtPart size

  size += 2 * part_size; // 2 * because we assumed 2 parts
  // End Ctxt content size

  // End eye-catcher
  size += 4;

  return size + offset;
}

int Squid::StorageOfOneElement(){
    if (!db_set){
        throw invalid_argument("ERROR: DB needs to be set to get storage cost");
    }
    return estimateCtxtSize(context, 0);
}
