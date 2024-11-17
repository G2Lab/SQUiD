#include <iostream>
#include <vector>
#include <fstream>

#include <benchmark/benchmark.h>

#include "../databases/FHE_SIMD_database.hpp"
#include "../databases/FHE_disk_database.hpp"

const int MOST_SNPS = 16;
const int MOST_SNPS_PRS = 16384;
const uint32_t SEED = 0; 
static FHEDiskDatabase *dbFHEInstance;

static void DoSetup(const benchmark::State &state)
{
    static bool callSetup = true;
    if (callSetup)
    {
        dbFHEInstance = new FHEDiskDatabase(constants::Large, "",true);
    }
    callSetup = false;
}

static void BM_CountQuery(benchmark::State &state)
{
    uint32_t commBytes = 0;
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    dbFHEInstance->setCompressedRows(state.range(2));

    bool conjunctive = state.range(1);
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }

    commBytes += query.size() * sizeof(query[0]);
    commBytes += sizeof(conjunctive);

    double start_capacity = dbFHEInstance->getDatabaseEntryCapacity();
    long start_level = dbFHEInstance->getDatabaseEntryLevel();
    
    double end_capacity;
    long end_level;

    for (auto _ : state)
    {
        auto result = dbFHEInstance->countQuery(conjunctive, query);
        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        end_capacity = result.bitCapacity();
        end_level = result.getPrimeSet().card();
        state.ResumeTiming();
        benchmark::DoNotOptimize(result);
    }
    commBytes += dbFHEInstance->storageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * dbFHEInstance->getSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
    state.counters["Start bit capacity"] = start_capacity;
    state.counters["End bit capacity"] = end_capacity;
    state.counters["Start level"] = start_level;
    state.counters["End level"] = end_level;
}

static void BM_MAFQuery(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(2));

    uint32_t commBytes = 0;
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    bool conjunctive = state.range(1);
    uint32_t snp = 0;

    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }

    commBytes += query.size() * sizeof(query[0]);
    commBytes += sizeof(conjunctive);
    commBytes += sizeof(snp);

    double start_capacity = dbFHEInstance->getDatabaseEntryCapacity();
    long start_level = dbFHEInstance->getDatabaseEntryLevel();
    
    double end_capacity;
    long end_level;

    for (auto _ : state)
    {
        auto result = dbFHEInstance->MAFQuery(snp, conjunctive, query);

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        end_capacity = result.bitCapacity();
        end_level = result.getPrimeSet().card();
        
        state.ResumeTiming();
        benchmark::DoNotOptimize(result);
    }
    commBytes += dbFHEInstance->storageOfOneElement();
    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * dbFHEInstance->getSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
    state.counters["Start bit capacity"] = start_capacity;
    state.counters["End bit capacity"] = end_capacity;
    state.counters["Start level"] = start_level;
    state.counters["End level"] = end_level;
}

static void BM_PRSQuery(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(1));

    uint32_t commBytes = 0;
    vector<pair<uint32_t, int32_t>> query = vector<pair<uint32_t, int32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(0, 0));
    }

    commBytes += query.size() * sizeof(query[0]);
    double start_capacity = dbFHEInstance->getDatabaseEntryCapacity();
    long start_level = dbFHEInstance->getDatabaseEntryLevel();

    double end_capacity;
    long end_level;

    for (auto _ : state)
    {
        auto result = dbFHEInstance->PRSQuery(query);
        state.PauseTiming();
        bool exceedError = false;
        for (size_t i = 0; i < result.size(); i++)
        {
            if (!result[i].isCorrect())
            {
                exceedError = true;
            }
        }
        if (exceedError)
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        end_capacity = result[0].bitCapacity();
        end_level = result[0].getPrimeSet().card();

        state.ResumeTiming();
        benchmark::DoNotOptimize(result);
    }
    commBytes += dbFHEInstance->getCompressedRows() * dbFHEInstance->storageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * dbFHEInstance->getSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
    state.counters["Start bit capacity"] = start_capacity;
    state.counters["End bit capacity"] = end_capacity;
    state.counters["Start level"] = start_level;
    state.counters["End level"] = end_level;
}

static void BM_SimilarityQuery(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(1));

    uint32_t commBytes = 0;

    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    for (int i = 0; i < state.range(0); i++)
    {
        d.push_back(dbFHEInstance->encrypt(0));
    }
    uint32_t targetSnp = 0;
    uint32_t threshold = 100;

    commBytes += d.size() * dbFHEInstance->storageOfOneElement();
    commBytes += sizeof(targetSnp);
    commBytes += sizeof(threshold);

    double start_capacity = dbFHEInstance->getDatabaseEntryCapacity();
    long start_level = dbFHEInstance->getDatabaseEntryLevel();

    double end_capacity;
    long end_level;

    for (auto _ : state)
    {
        auto result = dbFHEInstance->similarityQuery(targetSnp, d, threshold);

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        end_capacity = min(result.first.bitCapacity(), result.second.bitCapacity());
        end_level = result.first.getPrimeSet().card();

        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    commBytes += 2 * dbFHEInstance->storageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * dbFHEInstance->getSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
    state.counters["Start bit capacity"] = start_capacity;
    state.counters["End bit capacity"] = end_capacity;
    state.counters["Start level"] = start_level;
    state.counters["End level"] = end_level;
}


static void BM_CountQueryWithPKS(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(2));

    uint32_t commBytes = 0;
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    bool conjunctive = state.range(1);
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }

    const Meta& meta = dbFHEInstance->getMeta();
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(dbFHEInstance->getMeta().data->secretKey);
    std::vector<helib::DoubleCRT>& firstVectorRef = ksk.first;
    std::vector<helib::DoubleCRT>& secondVectorRef = ksk.second;

    std::pair<std::vector<helib::DoubleCRT>&, std::vector<helib::DoubleCRT>&> kskRef(firstVectorRef, secondVectorRef);

    commBytes += query.size() * sizeof(query[0]);
    commBytes += sizeof(conjunctive);

    for (auto _ : state)
    {
        auto result = dbFHEInstance->countQuery(conjunctive, query);
        result.PublicKeySwitch(kskRef);
        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();
        benchmark::DoNotOptimize(result);
    }
    commBytes += dbFHEInstance->storageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * dbFHEInstance->getSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
}

static void BM_MAFQueryWithPKS(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(2));

    uint32_t commBytes = 0;
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    bool conjunctive = state.range(1);
    uint32_t snp = 0;

    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }

    const Meta& meta = dbFHEInstance->getMeta();
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(dbFHEInstance->getMeta().data->secretKey);
    std::vector<helib::DoubleCRT>& firstVectorRef = ksk.first;
    std::vector<helib::DoubleCRT>& secondVectorRef = ksk.second;

    std::pair<std::vector<helib::DoubleCRT>&, std::vector<helib::DoubleCRT>&> kskRef(firstVectorRef, secondVectorRef);


    commBytes += query.size() * sizeof(query[0]);
    commBytes += sizeof(conjunctive);
    commBytes += sizeof(snp);

    for (auto _ : state)
    {
        auto result = dbFHEInstance->MAFQuery(snp, conjunctive, query);
        result.PublicKeySwitch(kskRef);

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    commBytes += dbFHEInstance->storageOfOneElement();
    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * dbFHEInstance->getSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
}

static void BM_PRSQueryWithPKS(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(1));
    uint32_t commBytes = 0;
    vector<pair<uint32_t, int32_t>> query = vector<pair<uint32_t, int32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(0, 0));
    }

    const Meta& meta = dbFHEInstance->getMeta();
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(dbFHEInstance->getMeta().data->secretKey);
    std::vector<helib::DoubleCRT>& firstVectorRef = ksk.first;
    std::vector<helib::DoubleCRT>& secondVectorRef = ksk.second;

    std::pair<std::vector<helib::DoubleCRT>&, std::vector<helib::DoubleCRT>&> kskRef(firstVectorRef, secondVectorRef);

    commBytes += query.size() * sizeof(query[0]);

    for (auto _ : state)
    {
        auto result = dbFHEInstance->PRSQuery(query);
        for (size_t i = 0; i < result.size(); i++)
        {
            result[i].PublicKeySwitch(kskRef);
        }
        state.PauseTiming();
        bool exceedError = false;
        for (size_t i = 0; i < result.size(); i++)
        {
            if (!result[i].isCorrect())
            {
                exceedError = true;
            }
        }
        if (exceedError)
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();
        benchmark::DoNotOptimize(result);
    }
    commBytes += dbFHEInstance->getCompressedRows() * dbFHEInstance->storageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * dbFHEInstance->getSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
}

static void BM_SimilarityQueryWithPKS(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(1));
    uint32_t commBytes = 0;

    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    for (int i = 0; i < state.range(0); i++)
    {
        d.push_back(dbFHEInstance->encrypt(0));
    }
    uint32_t targetSnp = 0;
    uint32_t threshold = 100;

    const Meta& meta = dbFHEInstance->getMeta();
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(dbFHEInstance->getMeta().data->secretKey);
    std::vector<helib::DoubleCRT>& firstVectorRef = ksk.first;
    std::vector<helib::DoubleCRT>& secondVectorRef = ksk.second;

    std::pair<std::vector<helib::DoubleCRT>&, std::vector<helib::DoubleCRT>&> kskRef(firstVectorRef, secondVectorRef);

    commBytes += d.size() * dbFHEInstance->storageOfOneElement();
    commBytes += sizeof(targetSnp);
    commBytes += sizeof(threshold);

    for (auto _ : state)
    {
        auto result = dbFHEInstance->similarityQuery(targetSnp, d, threshold);
        result.first.PublicKeySwitch(kskRef);
        result.second.PublicKeySwitch(kskRef);

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    commBytes += 2 * dbFHEInstance->storageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * dbFHEInstance->getSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
}


static void BM_RangeCountQuery(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(0));

    for (auto _ : state)
    {
        auto result = dbFHEInstance->countingRangeQuery(25, 75, 0);

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Number of patients"] = state.range(0) * dbFHEInstance->getSlotSize();
}

static void BM_RangeMAFQuery(benchmark::State &state)
{
    dbFHEInstance->setCompressedRows(state.range(0));

    for (auto _ : state)
    {
        auto result = dbFHEInstance->MAFRangeQuery(0, 25, 75, 0);

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Number of patients"] = state.range(0) * dbFHEInstance->getSlotSize();
}

BENCHMARK(BM_CountQuery)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, 1), {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_MAFQuery)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, 1), {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_PRSQuery)->ArgsProduct({{1024, 16384}, {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_SimilarityQuery)->ArgsProduct({{1024,16384}, {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_RangeCountQuery)->DenseRange(1, 10, 1)->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_RangeMAFQuery)->DenseRange(1, 10, 1)->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_CountQueryWithPKS)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, 1), {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_MAFQueryWithPKS)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, 1), {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_PRSQueryWithPKS)->ArgsProduct({{1024, 16384}, {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_SimilarityQueryWithPKS)->ArgsProduct({{1024, 16384}, {1,2,3,4,5,6,7,8,9,10}})->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK_MAIN();
