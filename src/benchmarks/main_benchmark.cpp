#include <iostream>
#include <vector>
#include <fstream>

#include <benchmark/benchmark.h>

#include "../server.hpp"

const int MOST_SNPS = 16;
const int MOST_SNPS_PRS = 16384;

static Server *serverInstance;

static void DoSetup(const benchmark::State &state)
{
    static bool callSetup = true;
    if (callSetup)
    {
        serverInstance = new Server(constants::BenchParams, true);

        serverInstance->GenData(1, MOST_SNPS);
    }
    callSetup = false;
}

static void BM_CountQuery(benchmark::State &state)
{
    if (state.range(2) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(2) - 1) * serverInstance->GetSlotSize(), MOST_SNPS);
    }

    uint32_t commBytes = 0;
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    bool conjunctive = state.range(1);
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }

    commBytes += query.size() * sizeof(query[0]);
    commBytes += sizeof(conjunctive);

    for (auto _ : state)
    {
        auto result = serverInstance->CountQuery(conjunctive, query);
        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();
        benchmark::DoNotOptimize(result);
    }
    commBytes += serverInstance->StorageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * serverInstance->GetSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
}

static void BM_MAFQuery(benchmark::State &state)
{
    if (state.range(2) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(2) - 1) * serverInstance->GetSlotSize(), MOST_SNPS);
    }
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

    for (auto _ : state)
    {
        auto result = serverInstance->MAFQuery(snp, conjunctive, query);

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    commBytes += serverInstance->StorageOfOneElement();
    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * serverInstance->GetSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
}

static void BM_PRSQuery(benchmark::State &state)
{
    if (state.range(1) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(1) - 1) * serverInstance->GetSlotSize(), 1);
    }
    uint32_t commBytes = 0;
    vector<pair<uint32_t, int32_t>> query = vector<pair<uint32_t, int32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(0, 0));
    }

    commBytes += query.size() * sizeof(query[0]);

    for (auto _ : state)
    {
        auto result = serverInstance->PRSQuery(query);
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
    commBytes += serverInstance->GetCompressedRows() * serverInstance->StorageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * serverInstance->GetSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
}

static void BM_SimilarityQuery(benchmark::State &state)
{
    if (state.range(1) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(1) - 1) * serverInstance->GetSlotSize(), MOST_SNPS);
    }
    uint32_t commBytes = 0;

    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    for (int i = 0; i < state.range(0); i++)
    {
        d.push_back(serverInstance->Encrypt(0));
    }
    uint32_t targetSnp = 0;
    uint32_t threshold = 100;

    commBytes += d.size() * serverInstance->StorageOfOneElement();
    commBytes += sizeof(targetSnp);
    commBytes += sizeof(threshold);

    for (auto _ : state)
    {
        auto result = serverInstance->SimilarityQuery(targetSnp, d, threshold);

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    commBytes += 2 * serverInstance->StorageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * serverInstance->GetSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
}

static void BM_CountQueryWithPKS(benchmark::State &state)
{
    if (state.range(2) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(2) - 1) * serverInstance->GetSlotSize(), MOST_SNPS);
    }

    uint32_t commBytes = 0;
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    bool conjunctive = state.range(1);
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }

    Meta meta;
    meta(constants::BenchParams);
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(serverInstance->GetMeta().data->secretKey);

    commBytes += query.size() * sizeof(query[0]);
    commBytes += sizeof(conjunctive);

    for (auto _ : state)
    {
        auto result = serverInstance->CountQuery(conjunctive, query);
        result.PublicKeySwitch(ksk);
        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();
        benchmark::DoNotOptimize(result);
    }
    commBytes += serverInstance->StorageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * serverInstance->GetSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
}

static void BM_MAFQueryWithPKS(benchmark::State &state)
{
    if (state.range(2) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(2) - 1) * serverInstance->GetSlotSize(), MOST_SNPS);
    }
    uint32_t commBytes = 0;
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    bool conjunctive = state.range(1);
    uint32_t snp = 0;

    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }

    Meta meta;
    meta(constants::BenchParams);
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(serverInstance->GetMeta().data->secretKey);

    commBytes += query.size() * sizeof(query[0]);
    commBytes += sizeof(conjunctive);
    commBytes += sizeof(snp);

    for (auto _ : state)
    {
        auto result = serverInstance->MAFQuery(snp, conjunctive, query);
        result.PublicKeySwitch(ksk);

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    commBytes += serverInstance->StorageOfOneElement();
    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(2) * serverInstance->GetSlotSize();
    state.counters["Number of filters"] = state.range(0);
    state.counters["Conjunctive (Or = 0, And = 1)"] = conjunctive;
}

static void BM_PRSQueryWithPKS(benchmark::State &state)
{
    if (state.range(1) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(1) - 1) * serverInstance->GetSlotSize(), 1);
    }
    uint32_t commBytes = 0;
    vector<pair<uint32_t, int32_t>> query = vector<pair<uint32_t, int32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(0, 0));
    }

    Meta meta;
    meta(constants::BenchParams);
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(serverInstance->GetMeta().data->secretKey);

    commBytes += query.size() * sizeof(query[0]);

    for (auto _ : state)
    {
        auto result = serverInstance->PRSQuery(query);
        for (size_t i = 0; i < result.size(); i++)
        {
            result[i].PublicKeySwitch(ksk);
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
    commBytes += serverInstance->GetCompressedRows() * serverInstance->StorageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * serverInstance->GetSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
}

static void BM_SimilarityQueryWithPKS(benchmark::State &state)
{
    if (state.range(1) != serverInstance->GetCompressedRows())
    {
        serverInstance->GenData(1 + (state.range(1) - 1) * serverInstance->GetSlotSize(), MOST_SNPS);
    }
    uint32_t commBytes = 0;

    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    for (int i = 0; i < state.range(0); i++)
    {
        d.push_back(serverInstance->Encrypt(0));
    }
    uint32_t targetSnp = 0;
    uint32_t threshold = 100;

    Meta meta;
    meta(constants::BenchParams);
    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);
    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(serverInstance->GetMeta().data->secretKey);

    commBytes += d.size() * serverInstance->StorageOfOneElement();
    commBytes += sizeof(targetSnp);
    commBytes += sizeof(threshold);

    for (auto _ : state)
    {
        auto result = serverInstance->SimilarityQuery(targetSnp, d, threshold);
        result.first.PublicKeySwitch(ksk);
        result.second.PublicKeySwitch(ksk);

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    commBytes += 2 * serverInstance->StorageOfOneElement();

    state.counters["Communication (B)"] = commBytes;
    state.counters["Number of patients"] = state.range(1) * serverInstance->GetSlotSize();
    state.counters["Number of SNPs"] = state.range(0);
}

static void BM_RangeCountQuery(benchmark::State &state)
{
    serverInstance->GenContinuousData(state.range(0) * serverInstance->GetSlotSize(), 1, 100);

    for (auto _ : state)
    {
        auto result = serverInstance->CountingRangeQuery(25, 75);

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Number of patients"] = state.range(0) * serverInstance->GetSlotSize();
}

static void BM_RangeMAFQuery(benchmark::State &state)
{
    serverInstance->GenContinuousData(state.range(0) * serverInstance->GetSlotSize(), 1, 100);
    serverInstance->GenData(state.range(0) * serverInstance->GetSlotSize(), 1);

    for (auto _ : state)
    {
        auto result = serverInstance->MAFRangeQuery(0, 25, 75);

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Number of patients"] = state.range(0) * serverInstance->GetSlotSize();
}

BENCHMARK(BM_RangeCountQuery)->DenseRange(1, 6, 1)->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_RangeMAFQuery)->DenseRange(1, 6, 1)->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_CountQuery)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, /*step=*/1), {1, 2, 3, 4, 5, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_MAFQuery)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, /*step=*/1), {1, 2, 3, 4, 5, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_PRSQuery)->ArgsProduct({{1024, 16384}, {1, 2, 3, 4, 5, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_SimilarityQuery)->ArgsProduct({{2, 16}, {1, 2, 3, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_CountQueryWithPKS)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, /*step=*/1), {1, 2, 3, 4, 5, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_MAFQueryWithPKS)->ArgsProduct({{2, 16}, benchmark::CreateDenseRange(0, 1, /*step=*/1), {1, 2, 3, 4, 5, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_PRSQueryWithPKS)->ArgsProduct({{1024, 16384}, {1, 2, 3, 4, 5, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_SimilarityQueryWithPKS)->ArgsProduct({{2, 16}, {1, 2, 3, 6}})->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK_MAIN();