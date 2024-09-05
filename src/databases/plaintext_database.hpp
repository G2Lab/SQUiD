#pragma once

#include <vector>
#include <string>
#include "database.hpp"
#include <random>
#include <iostream>

class PlaintextDatabase : public Database<int32_t, int32_t> {
public:
    PlaintextDatabase() = default;
    ~PlaintextDatabase() = default;

    void genData(uint32_t num_rows, uint32_t num_snp_cols, uint32_t seed) override;
    void genContinuousPhenoData(uint32_t num_cols, uint32_t low, uint32_t high, uint32_t seed) override;
    void genBinaryPhenoData(uint32_t num_cols, uint32_t seed) override;
    void setData(vector<vector<int32_t>> &db);
    void setBinaryPhenoData(vector<vector<int32_t>> &db);
    void setData(string vcf_file) override;

    int32_t getGenotype(uint32_t column, uint32_t row) const override;
    int32_t getContinuousPheno(uint32_t column, uint32_t row) const override;
    int32_t getBinaryPheno(uint32_t column, uint32_t row) const override;

    void setGenotype(int32_t data, uint32_t column, uint32_t compressed_row_index) override;

    int32_t countQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const override;
    int32_t MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const override;
    vector<int32_t> PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params) const override;
    pair<int32_t, int32_t> similarityQuery(uint32_t target_column, vector<int32_t> &d, uint32_t threshold) const override;

    int32_t countingRangeQuery(uint32_t lower, uint32_t upper, uint32_t phenotype_index) override;
    pair<int32_t, int32_t> MAFRangeQuery(uint32_t snp, uint32_t lower, uint32_t upper,uint32_t phenotype_index) override;

    void printData(bool with_headers) const override;
    void printMeta() const override;
};