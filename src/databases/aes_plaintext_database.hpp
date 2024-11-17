#pragma once

#include <vector>
#include <string>
#include "database.hpp"
#include <random>
#include <iostream>
#include <openssl/aes.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <openssl/err.h>
#include <cstring> 

class PlaintextAESDatabase : public Database<unsigned char*, int32_t> {
public:
    PlaintextAESDatabase();
    ~PlaintextAESDatabase();

    void initialise();

    void genData(uint32_t num_rows, uint32_t num_snp_cols, uint32_t seed) override;
    void genContinuousPhenoData(uint32_t num_cols, uint32_t low, uint32_t high, uint32_t seed) override;
    void genBinaryPhenoData(uint32_t num_cols, uint32_t seed) override;
    void setData(vector<vector<int32_t>> &db);
    void setBinaryPhenoData(vector<vector<int32_t>> &db);
    void setData(string vcf_file) override;

    void setGenotype(unsigned char* data, uint32_t column, uint32_t compressed_row_index);

    unsigned char* getGenotype(uint32_t column, uint32_t row) const;
    unsigned char* getContinuousPheno(uint32_t column, uint32_t row) const;
    unsigned char* getBinaryPheno(uint32_t column, uint32_t row) const;

    int32_t countQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const;
    int32_t MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const;
    vector<int32_t> PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params) const;
    pair<int32_t, int32_t> similarityQuery(uint32_t target_column, vector<unsigned char*> &d, uint32_t threshold) const;

    int32_t countingRangeQuery(uint32_t lower, uint32_t upper, uint32_t phenotype_index) override;
    pair<int32_t, int32_t> MAFRangeQuery(uint32_t snp, uint32_t lower, uint32_t upper,uint32_t phenotype_index) override;

    void printData(bool with_headers) const override;
    void printMeta() const override;

    unsigned char *key_128;
    unsigned char *iv;

    uint32_t num_slots;
};