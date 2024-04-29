#pragma once

#include <vector>
#include <helib/helib.h>
#include <string>
#include <sstream>
#include <map>


using namespace std;


class Squid
{
  public:
    Squid();

    void GenData(int _num_rows, int _num_cols);
    void SetData(vector<vector<unsigned long>> &db);
    void SetColumnHeaders(vector<string> &headers);

    void SetServerToExample();

    helib::Ctxt CountingQuery(bool conjunctive, vector<pair<int, int>>& query) const;
    pair<helib::Ctxt, helib::Ctxt> MAFQuery(int snp, bool conjunctive, vector<pair<int, int>> &query) const;
    vector<helib::Ctxt> PRSQuery(vector<pair<int, int>>& prs_params) const;
    vector<pair<helib::Ctxt, helib::Ctxt>> ChiSquareQuery(bool conjunctive, vector<pair<int, int>>& query, int disease_column, int number_of_chi);
    vector<pair<helib::Ctxt, helib::Ctxt>> ChiSquareQuery(int disease_column, int number_of_chi);

    vector<long> Decrypt(helib::Ctxt ctxt) const;

    void PrintContext() const;
    string PrintEncryptedDB(bool with_headers) const;
    int GetSlotSize() const;
    int GetNumRows() const;
    const helib::Context& GetContext() const;
    const helib::SecKey& GetSK() const;

    const vector<string>& GetColumnHeaders() const;


  private:
    void AddOneMod2(helib::Ctxt& a) const;
    helib::Ctxt MultiplyMany(vector<helib::Ctxt>& v) const;
    helib::Ctxt AddMany(vector<helib::Ctxt>& v) const;
    helib::Ctxt AddManySafe(vector<helib::Ctxt>& v) const;
    helib::Ctxt SquashCtxt(helib::Ctxt& ciphertext, int num_data_entries = 10) const;
    helib::Ctxt SquashCtxtLogTime(helib::Ctxt& ciphertext) const;
    helib::Ctxt SquashCtxtWithMask(helib::Ctxt& ciphertext, int index) const;
    helib::Ctxt EQTest(unsigned long a, const helib::Ctxt& b) const;
    vector<vector<helib::Ctxt>> filter(vector<pair<int, int>>& query) const;
    void CtxtExpand(helib::Ctxt &ciphertext) const;


    //Encrypt / Decrypt Methods
    helib::Ptxt<helib::BGV> DecryptPlaintext(helib::Ctxt ctxt) const;
    helib::Ctxt Encrypt(unsigned long a) const;
    helib::Ctxt Encrypt(vector<unsigned long> a) const;
    helib::Ctxt GetAnyElement() const;

    int StorageOfOneElement();
    const helib::Context context;
    helib::SecKey secret_key;
    helib::PubKey* public_key_ptr;
    
    bool db_set;
    
    int num_rows;
    int num_cols;
    int num_compressed_rows;
    int num_slots;
    
    vector<vector<helib::Ctxt>> encrypted_db; 
    std::map<std::string, std::pair<std::vector<helib::DoubleCRT>,std::vector<helib::DoubleCRT>>> key_switch_store;

    vector<string> column_headers;
    
    int one_over_two;
    int neg_three_over_two;
    int neg_one_over_two;

    int plaintext_modulus;
};

