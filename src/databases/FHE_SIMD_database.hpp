#pragma once

#include <vector>
#include <string>
#include "database.hpp"
#include <random>
#include <iostream>
#include <thread>
#include <utility>
#include <helib/helib.h>
#include "tools.hpp"
#include "comparator.hpp"
#include "../globals.hpp"

#define MAX_NUMBER_BITS 4
#define NOISE_THRES 0
#define EARLY_TERM_NOISE_THRES 50
#define WARN false

#define COMPRESSED 1

class FHESIMDDatabase : public Database<helib::Ctxt, helib::Ctxt>
{
public:
    FHESIMDDatabase(const Params &_params, bool _with_comparator);
    FHESIMDDatabase(const Params &_params, std::string context_filepath, bool _with_comparator);

    void initialize(const Params &_params, bool _with_comparator);

    ~FHESIMDDatabase() = default;

    void genData(uint32_t num_rows, uint32_t num_snp_cols, uint32_t seed) override;
    void genContinuousPhenoData(uint32_t num_continuous_pheno_cols, uint32_t low, uint32_t high, uint32_t seed) override;
    void genBinaryPhenoData(uint32_t num_binary_pheno_cols, uint32_t seed) override;

    void setData(vector<vector<int32_t>> &db) override;
    void setBinaryPhenoData(vector<vector<int32_t>> &db) override;
    void setData(string vcf_file) override;

    void setColumnHeaders(vector<string> &headers);

    helib::Ctxt getGenotype(uint32_t column, uint32_t row) const override;
    helib::Ctxt getContinuousPheno(uint32_t column, uint32_t row) const override;
    helib::Ctxt getBinaryPheno(uint32_t column, uint32_t row) const override;

    // Modify Operations
    void updateOneValue(uint32_t row, uint32_t col, uint32_t value);
    void updateOneRow(uint32_t row, vector<uint32_t> &vals);
    void insertOneRow(vector<uint32_t> &vals);
    void deleteRowAddition(uint32_t row);
    void deleteRowMultiplication(uint32_t row);

    void setGenotype(helib::Ctxt ctxt, uint32_t column, uint32_t compressed_row_index) override;

    // Querries
    helib::Ctxt countQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const override;
    helib::Ctxt countQueryP(vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads);

    helib::Ctxt MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const override;
    helib::Ctxt MAFQueryP(uint32_t snp, vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads);

    helib::Ctxt countingRangeQuery(uint32_t lower, uint32_t upper, uint32_t phenotype_index) override;
    pair<helib::Ctxt, helib::Ctxt> MAFRangeQuery(uint32_t snp, uint32_t lower, uint32_t upper, uint32_t phenotype_index) override;

    vector<helib::Ctxt> PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params) const override;
    helib::Ctxt PRSQueryP(vector<pair<uint32_t, int32_t>> &prs_params, uint32_t num_threads);

    pair<helib::Ctxt, helib::Ctxt> similarityQuery(uint32_t target_column, vector<helib::Ctxt> &d, uint32_t threshold) const override;
    pair<helib::Ctxt, helib::Ctxt> similarityQueryP(uint32_t target_column, vector<helib::Ctxt> &d, uint32_t num_threshold, uint32_t threads);

    helib::Ctxt countQueryPP(vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads);
    helib::Ctxt MAFQueryPP(uint32_t snp, vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads);
    vector<helib::Ctxt> PRSQueryPP(vector<pair<uint32_t, int32_t>> &prs_params, uint32_t num_threads);

    helib::Ctxt computeL2Norm(helib::Ctxt &d);
    helib::Ctxt computeOtherNorm(helib::Ctxt &d);

    helib::Ctxt squashCtxt(helib::Ctxt &ciphertext) const;
    helib::Ctxt squashCtxtLogTime(helib::Ctxt &ciphertext) const;
    helib::Ctxt squashCtxtLogTimePower2(helib::Ctxt &ciphertext) const;

    helib::Ctxt squashCtxtWithMask(helib::Ctxt &ciphertext, uint32_t index) const;
    void maskWithNumRows(vector<helib::Ctxt> &ciphertexts) const;
    helib::Ctxt EQTest(unsigned long a, helib::Ctxt b) const;
    vector<vector<helib::Ctxt>> filter(vector<pair<uint32_t, uint32_t>> &query) const;
    void ctxtExpand(helib::Ctxt &ciphertext) const;

    // Encrypt / Decrypt Methods
    helib::Ptxt<helib::BGV> decryptPlaintext(helib::Ctxt ctxt);
    vector<long> decrypt(helib::Ctxt ctxt) const;
    helib::Ctxt encrypt(unsigned long a);
    helib::Ctxt encryptFast(unsigned long a);
    helib::Ctxt encrypt(vector<unsigned long> a);
    helib::Ctxt encryptSK(unsigned long a);
    helib::Ctxt encryptSK(vector<unsigned long> a);
    helib::Ctxt getAnyElement();

    uint32_t getSlotSize() const;
    uint32_t getCompressedRows() const;
    uint32_t getCols() const;
    vector<string> getHeaders() const;
    const Meta &getMeta() const { return meta; }

    uint32_t storageOfOneElement() const;
    long getDatabaseEntryCapacity();
    long getDatabaseEntryLevel();

    int32_t getChunkingFactor() const;
    void setCompressedRows(uint32_t num_rows) { num_compressed_rows = num_rows; }
    void setNumberOfSNPColumns(uint32_t num_columns) { num_snp_cols = num_columns; }

    void printData(bool with_headers) const override;
    void printMeta() const override;

    unique_ptr<he_cmp::Comparator> comparator;

    vector<string> getColumnHeaders() const { return column_headers; }
    vector<string> getBinaryPhenoHeaders() const { return binary_pheno_headers; }
    vector<string> getContinuousPhenoHeaders() const { return continuous_pheno_headers; }

    void setBinaryPhenoHeaders(vector<string> headers) { binary_pheno_headers = headers; }
    void setContinuousPhenoHeaders(vector<string> headers) { continuous_pheno_headers = headers; }

    void printPhenoData(bool with_headers) const;

protected:
    Meta meta;

    uint32_t num_deletes = 0;
    uint32_t num_slots;

    bool with_comparator;

    uint32_t one_over_two;
    uint32_t neg_three_over_two;
    uint32_t neg_one_over_two;
    uint32_t neg_one;
    uint32_t neg_one_over_three;
    uint32_t neg_one_over_six;
    uint32_t one_over_six;

    uint32_t plaintext_modulus;

    helib::Ctxt *encryptZero;
    bool is_encrypt_zero_set = false;
    vector<string> binary_pheno_headers;
    vector<string> continuous_pheno_headers;
};