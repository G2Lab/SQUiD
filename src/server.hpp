#pragma once

#include <iostream>
#include <helib/helib.h>
#include <memory>
#include <string>
#include <vector>
#include "globals.hpp"
#include "comparator.hpp"
#include "tools.hpp"
#include <thread>
#include <utility>

#define MAX_NUMBER_BITS 4
#define NOISE_THRES 2
#define WARN false

using namespace std;


class Server{
public:
    
    //Setup
    Server(const Params& _params, bool _using_disk);
    
    void GenData(uint32_t  _num_rows, uint32_t  _num_cols);  
    void GenContinuousData(uint32_t _num_rows, uint32_t _low, uint32_t _high);
    void GenDataDummy(uint32_t  _num_rows, uint32_t  _num_cols);    
    void SetData(vector<vector<uint32_t >> &db);    
    void SetData(string vcf_file);

    void SetColumnHeaders(vector<string> &headers);
    
    //Modify Operations
    void UpdateOneValue(uint32_t  row, uint32_t  col, uint32_t  value);
    void UpdateOneRow(uint32_t  row, vector<uint32_t > &vals);
    void InsertOneRow(vector<uint32_t > &vals);
    void DeleteRowAddition(uint32_t  row);
    void DeleteRowMultiplication(uint32_t  row);
    
    //Querries
    helib::Ctxt CountQuery(bool conjunctive, vector<pair<uint32_t , uint32_t >>& query);
    helib::Ctxt CountQueryP(vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads);
    helib::Ctxt MAFQuery(uint32_t  snp, bool conjunctive, vector<pair<uint32_t , uint32_t >> &query);
    helib::Ctxt MAFQueryP(uint32_t  snp, vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads);

    helib::Ctxt CountingRangeQuery(uint32_t  lower, uint32_t  upper);
    pair<helib::Ctxt, helib::Ctxt> MAFRangeQuery(uint32_t  snp, uint32_t  lower, uint32_t  upper);

    vector<helib::Ctxt> PRSQuery(vector<pair<uint32_t , int32_t >>& prs_params);
    helib::Ctxt PRSQueryP(vector<pair<uint32_t, int32_t>> &prs_params, uint32_t num_threads);
    pair<helib::Ctxt, helib::Ctxt> SimilarityQuery(uint32_t  target_column, vector<helib::Ctxt>& d, uint32_t  threshold);
    pair<helib::Ctxt, helib::Ctxt> SimilarityQueryP(uint32_t  target_column, vector<helib::Ctxt>& d, uint32_t  num_threshold, uint32_t  threads);

    void AddOneMod2(helib::Ctxt& a);
    helib::Ctxt SquashCtxt(helib::Ctxt& ciphertext, uint32_t  num_data_entries = 10);
    helib::Ctxt SquashCtxtLogTime(helib::Ctxt& ciphertext);
    helib::Ctxt SquashCtxtLogTimePower2(helib::Ctxt& ciphertext);

    helib::Ctxt SquashCtxtWithMask(helib::Ctxt& ciphertext, uint32_t  index);
    void MaskWithNumRows(vector<helib::Ctxt>& ciphertexts);
    helib::Ctxt EQTest(unsigned long a, helib::Ctxt& b);
    vector<vector<helib::Ctxt>> filter(vector<pair<uint32_t , uint32_t >>& query);
    void CtxtExpand(helib::Ctxt &ciphertext);
    
    //Encrypt / Decrypt Methods
    helib::Ptxt<helib::BGV> DecryptPlaintext(helib::Ctxt ctxt);
    vector<long> Decrypt(helib::Ctxt ctxt);
    helib::Ctxt Encrypt(unsigned long a);
    helib::Ctxt Encrypt(vector<unsigned long> a);
    helib::Ctxt EncryptSK(unsigned long a);
    helib::Ctxt EncryptSK(vector<unsigned long> a);
    helib::Ctxt GetAnyElement();
    
    void PrintContext();
    void PrintEncryptedDB(bool with_headers);
    uint32_t GetSlotSize();
    uint32_t GetCompressedRows();
    uint32_t GetCols();
    vector<string> GetHeaders();
    Meta& GetMeta(){return meta;}

    uint32_t  StorageOfOneElement();
    
private:
    Meta meta;

    unique_ptr<he_cmp::Comparator> comparator;
    
    bool db_set;
    bool with_similarity;
    
    uint32_t  num_rows = 0;
    uint32_t  num_cols = 0;
    uint32_t  num_compressed_rows = 0;
    uint32_t  num_slots;
    uint32_t  num_deletes = 0;
    
    vector<vector<helib::Ctxt>> encrypted_db; 
    vector<string> column_headers;

    vector<helib::Ctxt> continuous_db;
    
    uint32_t  one_over_two;
    uint32_t  neg_three_over_two;
    uint32_t  neg_one_over_two;

    uint32_t  plaintext_modulus;
};