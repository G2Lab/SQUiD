#include "FHE_SIMD_database.hpp"

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------

//                                                         HELPER FUNCTIONS

// ------------------------------------------------------------------------------------------------------------------------

helib::Ctxt AddMany(vector<helib::Ctxt> &v);
helib::Ctxt AddManySafe(vector<helib::Ctxt> &v, const PubKey &pk);
helib::Ctxt MultiplyMany(vector<helib::Ctxt> &v);
void AddOneMod2(helib::Ctxt &a);

FHESIMDDatabase::FHESIMDDatabase(const Params &_params, bool _with_similarity)
{
    meta(_params);
    initialize(_params, _with_similarity);
}

FHESIMDDatabase::FHESIMDDatabase(const Params &_params, string context_filepath, bool _with_similarity)
{
    meta(context_filepath, _params);
    initialize(_params, _with_similarity);
}

void FHESIMDDatabase::initialize(const Params &_params, bool _with_similarity)
{
    num_slots = meta.data->ea.size();
    plaintext_modulus = meta.data->context.getP();

    one_over_two = get_inverse(1, 2, plaintext_modulus);
    neg_three_over_two = get_inverse(-3, 2, plaintext_modulus);
    neg_one_over_two = get_inverse(-1, 2, plaintext_modulus);

    neg_one = plaintext_modulus - 1;
    neg_one_over_three = get_inverse(-1, 3, plaintext_modulus);
    neg_one_over_six = get_inverse(-1, 6, plaintext_modulus);
    one_over_six = get_inverse(1, 6, plaintext_modulus);

    snp_data_set = false;

    with_comparator = _with_similarity;

    if (with_comparator)
    {
        comparator = std::make_unique<he_cmp::Comparator>(meta.data->context, he_cmp::UNI, _params.d, _params.l, meta.data->secretKey, false);
    }

    encryptZero = new helib::Ctxt(meta.data->publicKey);
}

helib::Ctxt FHESIMDDatabase::getGenotype(uint32_t column, uint32_t row) const
{
    if (column >= num_snp_cols || row >= num_compressed_rows || !snp_data_set)
    {
        throw invalid_argument("ERROR: Trying to access out of bounds SNP data or SNP data not set");
    }
    return snp_data[column][row];
}

helib::Ctxt FHESIMDDatabase::getContinuousPheno(uint32_t column, uint32_t row) const
{
    if (column >= num_continuous_pheno_cols || row >= num_compressed_rows || !continuous_phenotype_data_set)
    {
        throw invalid_argument("ERROR: Trying to access out of bounds continuous phenotype data or continuous phenotype data not set");
    }
    return continuous_phenotype_data[column][row];
}

helib::Ctxt FHESIMDDatabase::getBinaryPheno(uint32_t column, uint32_t row) const
{
    if (column >= num_binary_pheno_cols || row >= num_compressed_rows || !binary_phenotype_data_set)
    {
        throw invalid_argument("ERROR: Trying to access out of bounds binary phenotype data or binary phenotype data not set");
    }
    return binary_phenotype_data[column][row];
}

helib::Ctxt FHESIMDDatabase::EQTest(unsigned long a, helib::Ctxt b) const
{
    helib::Ctxt clone = b;
    DynamicCtxtPowers babyStep(clone, 3);

    switch (a)
    {
    case 0:
    {
        // f(x) = x^3 / 2 - x^2 - x / 2 + 1
        // -1 -> 0
        //  0 -> 1
        //  1 -> 0
        //  2 -> 0

        NTL::ZZX poly;
        NTL::SetCoeff(poly, 3, one_over_two);
        NTL::SetCoeff(poly, 2, neg_one);
        NTL::SetCoeff(poly, 1, neg_one_over_two);
        NTL::SetCoeff(poly, 0, 1);

        simplePolyEval(clone, poly, babyStep);

        return clone;
    }
    case 1:
    {
        // f(x) = -x^3 / 2 + x^2 / 2 + x
        // -1 -> 0
        //  0 -> 0
        //  1 -> 1
        //  2 -> 0
        NTL::ZZX poly;
        NTL::SetCoeff(poly, 3, neg_one_over_two);
        NTL::SetCoeff(poly, 2, one_over_two);
        NTL::SetCoeff(poly, 1, 1);

        simplePolyEval(clone, poly, babyStep);

        return clone;
    }
    case 2:
    {
        // f(x) = x^3 / 6 - x / 6
        // -1 -> 0
        //  0 -> 0
        //  1 -> 0
        //  2 -> 1

        NTL::ZZX poly;
        NTL::SetCoeff(poly, 3, one_over_six);
        NTL::SetCoeff(poly, 1, neg_one_over_six);

        simplePolyEval(clone, poly, babyStep);

        return clone;
    }
    default:
        cout << "Can't use a value of a other than 0, 1, or 2" << endl;
        throw invalid_argument("ERROR: invalid value for EQTest");
    }
}

vector<vector<helib::Ctxt>> FHESIMDDatabase::filter(vector<pair<uint32_t, uint32_t>> &query) const
{
    vector<vector<helib::Ctxt>> feature_cols;

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        vector<helib::Ctxt> indv_vector;
        for (pair<uint32_t, uint32_t> i : query)
        {
            indv_vector.push_back(EQTest(i.second, getGenotype(i.first, j)));

            if (constants::DEBUG == 3)
            {
                cout << "checking equality to " << i.second << endl;
                cout << "original:";
                print_vector(decrypt(getGenotype(i.first, j)));
                cout << "result  :";
                print_vector(decrypt(EQTest(i.second, getGenotype(i.first, j))));
            }
        }
        feature_cols.push_back(indv_vector);
    }
    return feature_cols;
}

vector<long> FHESIMDDatabase::decrypt(helib::Ctxt ctxt) const
{
    if (ctxt.capacity() < NOISE_THRES)
    {
        cout << "NOISE BOUNDS EXCEEDED!!!" << endl;
    }
    helib::Ptxt<helib::BGV> new_plaintext_result(meta.data->context);
    meta.data->secretKey.Decrypt(new_plaintext_result, ctxt);

    vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();

    vector<long> result = vector<long>(num_slots);

    for (uint32_t i = 0; i < num_slots; i++)
    {
        result[i] = (long)poly_mod_result[i];
    }

    return result;
}

helib::Ptxt<helib::BGV> FHESIMDDatabase::decryptPlaintext(helib::Ctxt ctxt)
{

    if (constants::DEBUG && ctxt.capacity() < NOISE_THRES)
    {
        cout << "NOISE BOUNDS EXCEEDED!!!" << endl;
    }

    helib::Ptxt<helib::BGV> new_plaintext_result(meta.data->context);
    meta.data->secretKey.Decrypt(new_plaintext_result, ctxt);

    return new_plaintext_result;
}

helib::Ctxt FHESIMDDatabase::encrypt(unsigned long a)
{
    vector<unsigned long> a_vec = vector<unsigned long>();
    for (size_t i = 0; i < num_slots; i++)
    {
        a_vec.push_back(a);
    }

    return encrypt(a_vec);
}

helib::Ctxt FHESIMDDatabase::encryptFast(unsigned long a)
{
    helib::Ptxt<helib::BGV> ptxt(meta.data->context);
    for (size_t i = 0; i < num_slots; i++)
    {
        ptxt[i] = a;
    }
    helib::Ctxt ctxt = *encryptZero;
    ctxt += ptxt;
    return ctxt;
}

helib::Ctxt FHESIMDDatabase::encrypt(vector<unsigned long> a)
{
    if (a.size() > num_slots)
    {
        throw invalid_argument("Trying to encrypt vector with too many elements");
    }
    helib::Ptxt<helib::BGV> ptxt(meta.data->context);

    for (size_t i = 0; i < a.size(); ++i)
    {
        ptxt[i] = a[i];
    }

    helib::Ctxt ctxt(meta.data->publicKey);

    meta.data->publicKey.Encrypt(ctxt, ptxt);

    return ctxt;
}

helib::Ctxt FHESIMDDatabase::encryptSK(unsigned long a)
{
    vector<unsigned long> a_vec = vector<unsigned long>();
    for (size_t i = 0; i < num_slots; i++)
    {
        a_vec.push_back(a);
    }

    return encryptSK(a_vec);
}

helib::Ctxt FHESIMDDatabase::encryptSK(vector<unsigned long> a)
{
    if (a.size() > num_slots)
    {
        throw invalid_argument("Trying to encrypt vector with too many elements");
    }
    helib::Ptxt<helib::BGV> ptxt(meta.data->context);

    for (size_t i = 0; i < a.size(); ++i)
    {
        ptxt[i] = a[i];
    }

    helib::Ctxt ctxt(meta.data->secretKey);

    EncodedPtxt eptxt;
    ptxt.encode(eptxt);

    meta.data->secretKey.Encrypt(ctxt, eptxt);

    return ctxt;
}

helib::Ctxt FHESIMDDatabase::getAnyElement()
{
    return getGenotype(0, 0);
}

void FHESIMDDatabase::printMeta() const
{
    meta.data->context.printout();
    cout << endl;
    cout << "Security: " << meta.data->context.securityLevel() << endl;
    cout << "Num slots: " << num_slots << endl;
    cout << "Num rows:" << num_rows << endl;
    cout << "Num compressed rows: " << num_compressed_rows << endl;

    cout << "Number of snp columns: " << num_snp_cols << endl;
    cout << "Number of binary pheno columns: " << num_binary_pheno_cols << endl;
    cout << "Number of continuous pheno columns: " << num_continuous_pheno_cols << endl;
}

void FHESIMDDatabase::printData(bool with_headers) const
{
    if (with_headers)
    {
        vector<uint32_t> string_length_count = vector<uint32_t>();

        cout << "|";
        for (uint32_t i = 0; i < num_snp_cols; i++)
        {
            cout << column_headers[i] << "|";
            string_length_count.push_back(column_headers[i].length());
        }
        cout << endl;
        cout << "--------------";
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {

            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (uint32_t i = 0; i < num_snp_cols; i++)
            {
                temp_storage.push_back(decrypt(getGenotype(i, j)));
            }
            for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
            {

                cout << endl
                     << "|";

                for (uint32_t i = 0; i < num_snp_cols; i++)
                {
                    if (i > 0)
                    {
                        for (uint32_t space = 0; space < string_length_count[i]; space++)
                        {
                            cout << " ";
                        }
                    }

                    cout << temp_storage[i][jj];

                    if (i == num_snp_cols - 1)
                    {
                        for (uint32_t space = 0; space < string_length_count[i]; space++)
                        {
                            cout << " ";
                        }
                        cout << "|";
                    }
                }
            }
        }
    }
    else
    {
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {

            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (uint32_t i = 0; i < num_snp_cols; i++)
            {
                temp_storage.push_back(decrypt(getGenotype(i, j)));
            }
            for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
            {

                cout << endl
                     << "|";

                for (uint32_t i = 0; i < num_snp_cols; i++)
                {
                    if (i > 0)
                    {
                        cout << " ";
                    }

                    cout << temp_storage[i][jj];

                    if (i == num_snp_cols - 1)
                    {
                        cout << "|";
                    }
                }
            }
        }
    }
    cout << endl;
}

void FHESIMDDatabase::printPhenoData(bool with_headers) const
{
    if (binary_phenotype_data_set)
    {
        if (with_headers)
        {
            vector<uint32_t> string_length_count = vector<uint32_t>();

            std::cout << "|";
            for (uint32_t i = 0; i < num_binary_pheno_cols; i++)
            {
                std::cout << binary_pheno_headers[i] << "|";
                string_length_count.push_back(binary_pheno_headers[i].length());
            }
            std::cout << std::endl;
            std::cout << "--------------";
            for (uint32_t j = 0; j < num_compressed_rows; j++)
            {
                vector<vector<long>> temp_storage = vector<vector<long>>();
                for (uint32_t i = 0; i < num_binary_pheno_cols; i++)
                {
                    temp_storage.push_back(decrypt(getBinaryPheno(i, j)));
                }
                for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
                {

                    std::cout << std::endl
                              << "|";

                    for (uint32_t i = 0; i < num_binary_pheno_cols; i++)
                    {
                        if (i > 0)
                        {
                            for (uint32_t space = 0; space < string_length_count[i]; space++)
                            {
                                std::cout << " ";
                            }
                        }

                        std::cout << temp_storage[i][jj];

                        if (i == num_binary_pheno_cols - 1)
                        {
                            for (uint32_t space = 0; space < string_length_count[i]; space++)
                            {
                                std::cout << " ";
                            }
                            std::cout << "|";
                        }
                    }
                }
            }
        }
        else
        {
            throw runtime_error("Not implemented");
        }
    }
    std::cout << std::endl;
    if (continuous_phenotype_data_set)
    {
        if (with_headers)
        {
            vector<uint32_t> string_length_count = vector<uint32_t>();

            std::cout << "|";
            for (uint32_t i = 0; i < num_continuous_pheno_cols; i++)
            {
                std::cout << continuous_pheno_headers[i] << "|";
                string_length_count.push_back(continuous_pheno_headers[i].length());
            }
            std::cout << std::endl;
            std::cout << "--------------";
            for (uint32_t j = 0; j < num_compressed_rows; j++)
            {
                vector<vector<long>> temp_storage = vector<vector<long>>();
                for (uint32_t i = 0; i < num_continuous_pheno_cols; i++)
                {
                    temp_storage.push_back(decrypt(getContinuousPheno(i, j)));
                }
                for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
                {

                    std::cout << std::endl
                              << "|";

                    for (uint32_t i = 0; i < num_continuous_pheno_cols; i++)
                    {
                        if (i > 0)
                        {
                            for (uint32_t space = 0; space < string_length_count[i]; space++)
                            {
                                std::cout << " ";
                            }
                        }

                        std::cout << temp_storage[i][jj];

                        if (i == num_continuous_pheno_cols - 1)
                        {
                            for (uint32_t space = 0; space < string_length_count[i]; space++)
                            {
                                std::cout << " ";
                            }
                            std::cout << "|";
                        }
                    }
                }
            }
            std::cout << std::endl;
        }
        else
        {
            throw runtime_error("Not implemented");
        }
    }
}

vector<string> FHESIMDDatabase::getHeaders() const
{
    return column_headers;
}

// IMPORTED FROM HELIB SOURCE CODE
inline long estimateCtxtSize(const helib::Context &context, long offset)
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

uint32_t FHESIMDDatabase::storageOfOneElement() const
{
    if (!snp_data_set)
    {
        throw invalid_argument("ERROR: DB needs to be set to get storage cost");
    }

    return estimateCtxtSize(meta.data->context, 0);
}

long FHESIMDDatabase::getDatabaseEntryCapacity()
{
    if (!snp_data_set)
    {
        throw invalid_argument("ERROR: DB needs to be set to get storage cost");
    }
    return getAnyElement().bitCapacity();
}

long FHESIMDDatabase::getDatabaseEntryLevel()
{
    if (!snp_data_set)
    {
        throw invalid_argument("ERROR: DB needs to be set to get storage cost");
    }
    return getAnyElement().getPrimeSet().card();
}

uint32_t FHESIMDDatabase::getCols() const
{
    return num_snp_cols;
}
uint32_t FHESIMDDatabase::getSlotSize() const
{
    return num_slots;
}
uint32_t FHESIMDDatabase::getCompressedRows() const
{
    return num_compressed_rows;
}

int32_t FHESIMDDatabase::getChunkingFactor() const
{
    int32_t depth = floor(log2(min(num_slots, num_rows)));
    uint32_t p = meta.data->params.p;

    for (int32_t i = 0; i < depth; i++)
    {
        uint32_t max_size = (long)pow(2, i + 2) * num_compressed_rows;
        if (max_size > p)
        {
            return i;
        }
    }
    return -1;
}

template <typename T, typename Allocator>
void print_vector(const vector<T, Allocator> &vect, int num_entries)
{
    cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++)
    {
        cout << ", " << vect[i];
    }
    cout << endl;
}
