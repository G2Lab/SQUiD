#pragma once

#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <thread>
#include <utility>
#include <list>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <memory>
#include <helib/helib.h>
#include "FHE_disk_database.hpp"
#include "FHE_SIMD_database.hpp"
#include "tools.hpp"
#include "comparator.hpp"
#include "../globals.hpp"

#define DISK_DIR "/data/db_disk"
#define META_FILEPATH "meta"
#define DB_META_FILE "meta.db"

class FHEDiskDatabase : public FHESIMDDatabase
{
public:
    // starting up disk database for first time
    FHEDiskDatabase(const Params &_params, bool _with_similarity)
        : FHESIMDDatabase(_params, _with_similarity)
    {
        DISK_DIR_FULL = getProjectRootPath() + DISK_DIR;
        if (!checkIfDiskDirExists())
        {
            createDiskDir();
        }
        storeMetadata();
    }
    // starting up disk database for subsequent runs
    FHEDiskDatabase(const Params &_params, std::string context_filepath, bool _with_similarity) : FHESIMDDatabase(_params, getProjectRootPath() + DISK_DIR + "/" + META_FILEPATH, _with_similarity)
    {
        DISK_DIR_FULL = getProjectRootPath() + DISK_DIR;

        if (!checkIfDiskDirExists())
        {
            createDiskDir();
            storeMetadata();
        }
        loadDBMetadata();
    }
    ~FHEDiskDatabase() = default;

    void genData(uint32_t num_rows, uint32_t num_snp_cols, uint32_t seed) override;
    void genBinaryPhenoData(uint32_t num_binary_pheno_cols, uint32_t seed) override;
    void genContinuousPhenoData(uint32_t num_continuous_pheno_cols, uint32_t low, uint32_t high, uint32_t seed) override;

    void setVCFData(string vcf_file);
    void setPhenoData(string pheno_file);

    void setData(vector<vector<int32_t>> &db) override;
    void multithreadSetData(string vcf_file, uint32_t num_threads);
    void setNumThreadsEncryption(uint32_t threads) { num_threads_encryption = threads; }

    void setGenotype(helib::Ctxt ctxt, uint32_t column, uint32_t compressed_row_index) override;

    helib::Ctxt getGenotype(uint32_t column, uint32_t row) const override;
    helib::Ctxt getContinuousPheno(uint32_t column, uint32_t row) const override;
    helib::Ctxt getBinaryPheno(uint32_t column, uint32_t row) const override;

    void storeMetadata();

    void storeDBMetadata();
    void loadDBMetadata();

    bool checkIfDiskDirExists();
    void createDiskDir();

    void createColumnDir(uint32_t column_idx);
    void createColumnDirPheno(uint32_t column_idx, std::string postfix);
    void saveCtxt(helib::Ctxt &ctxt, uint32_t column_idx, uint32_t row_idx);
    void saveCtxtPheno(helib::Ctxt &ctxt, uint32_t column_idx, uint32_t row_idx, std::string postfix);

    std::string DISK_DIR_FULL;

    // API Helping functions
    string printDBString(bool with_headers) const;
    const helib::Context &getContext() const { return meta.data->context; }
    const helib::SecKey &getSecretKey() const { return meta.data->secretKey; }
    uint32_t getNumRows() const { return num_rows; }
private:
    // Define the cache structures
    uint32_t num_threads_encryption = 1;

};