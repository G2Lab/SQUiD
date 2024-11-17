#pragma once

#include <vector>
#include <string>

using namespace std;
template <typename T, typename R>
class Database
{
public:
    virtual ~Database() = default;

    virtual void genData(uint32_t num_rows, uint32_t num_snp_cols, uint32_t seed) = 0;
    virtual void genContinuousPhenoData(uint32_t num_cols, uint32_t low, uint32_t high, uint32_t seed) = 0;
    virtual void genBinaryPhenoData(uint32_t num_cols, uint32_t seed) = 0;
    virtual void setData(vector<vector<int32_t>> &db) = 0;
    virtual void setBinaryPhenoData(vector<vector<int32_t>> &db) = 0;
    virtual void setData(string vcf_file) = 0;

    virtual T getGenotype(uint32_t column, uint32_t row) const = 0;
    virtual T getContinuousPheno(uint32_t column, uint32_t row) const = 0;
    virtual T getBinaryPheno(uint32_t column, uint32_t row) const = 0;

    virtual void setGenotype(T data, uint32_t column, uint32_t compressed_row_index) = 0;

    virtual R countQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const = 0;
    virtual R MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const = 0;
    virtual vector<R> PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params) const = 0;
    virtual pair<R, R> similarityQuery(uint32_t target_column, vector<T> &d, uint32_t threshold) const = 0;

    virtual R countingRangeQuery(uint32_t lower, uint32_t upper, uint32_t phenotype_index) = 0;
    virtual pair<R, R> MAFRangeQuery(uint32_t snp, uint32_t lower, uint32_t upper, uint32_t phenotype_index) = 0;

    uint32_t num_snp_cols = 0;
    uint32_t num_continuous_pheno_cols = 0;
    uint32_t num_binary_pheno_cols = 0;

    uint32_t num_rows;
    uint32_t num_compressed_rows;

    virtual void printData(bool with_headers) const = 0;
    virtual void printMeta() const = 0;

    bool snp_data_set;
    bool continuous_phenotype_data_set;
    bool binary_phenotype_data_set;
    bool headers_set;

    vector<vector<T>> snp_data;
    vector<string> column_headers;
    vector<vector<T>> continuous_phenotype_data;
    vector<vector<T>> binary_phenotype_data;
};