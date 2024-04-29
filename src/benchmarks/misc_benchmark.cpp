#include <iostream>
#include <vector>
#include <fstream>

#include <benchmark/benchmark.h>

#include "../server.hpp"

static Server *serverInstance;

static const unsigned long Q = 65536;
static const unsigned long P = 131071;
static const unsigned long R = 1;
static const unsigned long BITS = 880;

static void DoSetup(const benchmark::State &state)
{
    static bool callSetup = true;
    if (callSetup)
    {
        serverInstance = new Server(constants::BenchParams, true);
        serverInstance->PrintContext();
    }
    serverInstance->GenData(1, 16);
    callSetup = false;
}

static void BM_EncrpytCiphertext(benchmark::State &state)
{
    std::vector<unsigned long> data = std::vector<unsigned long>(1000);

    for (auto _ : state)
    {
        auto result = serverInstance->Encrypt(data);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_StorageCiphertext(benchmark::State &state)
{
    std::vector<unsigned long> data = std::vector<unsigned long>(1);

    uint32_t result = 0;

    for (auto _ : state)
    {
        auto ctxt = serverInstance->Encrypt(data);
        std::ofstream outFile("tempfile.bin", std::ios::binary);
        ctxt.writeTo(outFile);
        outFile.close();

        std::ifstream inFile("tempfile.bin", std::ios::binary);

        inFile.seekg(0, std::ios::end);
        std::streampos fileSize = inFile.tellg();
        inFile.close();

        result = static_cast<int>(fileSize);

        benchmark::DoNotOptimize(ctxt);
        benchmark::DoNotOptimize(result);
    }
    state.counters["Storage (B)"] = result;
    state.counters["Number of patients"] = 1;
}

static void BM_GeneratePublicKeySwitch(benchmark::State &state)
{
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                                 .m(Q)
                                 .p(P)
                                 .r(R)
                                 .bits(BITS)
                                 .c(state.range(0))
                                 .build();
    uint32_t size = 0;

    helib::SecKey owner_secret_key(context);
    owner_secret_key.GenSecKey();
    helib::PubKey owner_public_key(owner_secret_key);

    helib::SecKey client_secret_key(context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);

    for (auto _ : state)
    {
        pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(owner_secret_key, 0);
        benchmark::DoNotOptimize(ksk);

        state.PauseTiming();
        std::stringstream stream;
        ksk.first.at(0).writeTo(stream);
        std::string serializedData = stream.str();
        uint32_t numBytes = serializedData.size();
        size = (ksk.first.size() + ksk.second.size()) * numBytes;
        state.ResumeTiming();
    }

    state.counters["Storage (B)"] = size;
}

static void BM_SwithPublicKeySwitch(benchmark::State &state)
{
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                                 .m(Q)
                                 .p(P)
                                 .r(R)
                                 .bits(BITS)
                                 .c(state.range(0))
                                 .build();

    helib::SecKey owner_secret_key(context);
    owner_secret_key.GenSecKey();
    helib::PubKey owner_public_key(owner_secret_key);

    helib::SecKey client_secret_key(context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);

    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(owner_secret_key, 0);

    helib::Ptxt<helib::BGV> ptxt(context);
    int num_slots = 10;
    vector<long> original_values = vector<long>(num_slots);
    for (int i = 0; i < num_slots; i++)
    {
        ptxt[i] = i;
        original_values[i] = i;
    }

    helib::Ctxt ctxt(owner_public_key);
    owner_public_key.Encrypt(ctxt, ptxt);

    for (auto _ : state)
    {
        state.PauseTiming();
        helib::Ctxt clone = ctxt;
        state.ResumeTiming();

        clone.PublicKeySwitch(ksk.first, ksk.second, 0, client_secret_key);
        benchmark::DoNotOptimize(clone);
    }
}

static void BM_ParallelSimilarityQuery(benchmark::State &state)
{
    helib::Ctxt d = serverInstance->Encrypt(0);
    uint32_t targetSnp = 0;
    uint32_t threshold = 100;

    for (auto _ : state)
    {
        auto result = serverInstance->SimilarityQueryP(targetSnp, d, threshold, state.range(1), state.range(0));

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
}

static void BM_ParallelPRSQuery(benchmark::State &state)
{
    vector<pair<uint32_t, int32_t>> query = vector<pair<uint32_t, int32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(0, 0));
    }
    for (auto _ : state)
    {
        auto result = serverInstance->PRSQueryP(query, state.range(1));

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
}

static void BM_ParallelMAFQuery(benchmark::State &state)
{
    if (state.range(0) > serverInstance->GetCols())
    {
        serverInstance->GenData(1, state.range(0));
    }

    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }
    for (auto _ : state)
    {
        auto result = serverInstance->MAFQueryP(0, query, state.range(1));

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
}

static void BM_ParallelCountQuery(benchmark::State &state)
{
    if (state.range(0) > serverInstance->GetCols())
    {
        serverInstance->GenData(1, state.range(0));
    }

    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }
    for (auto _ : state)
    {
        auto result = serverInstance->CountQueryP(query, state.range(1));

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
}

static void BM_UpdateOneValue(benchmark::State &state)
{
    int db_snps = state.range(0);
    serverInstance->GenData(1, db_snps);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        serverInstance->UpdateOneValue(0, 0, 0);
    }
}
static void BM_UpdateOneRow(benchmark::State &state)
{
    int db_snps = state.range(0);
    serverInstance->GenData(1, db_snps);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        serverInstance->UpdateOneRow(0, vals);
    }
}

static void BM_InsertRow(benchmark::State &state)
{
    int db_snps = state.range(0);
    serverInstance->GenData(1, db_snps);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        serverInstance->InsertOneRow(vals);
    }
}
static void BM_DeleteRowAddition(benchmark::State &state)
{
    int db_snps = state.range(0);
    serverInstance->GenData(1, db_snps);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        serverInstance->DeleteRowAddition(0);
    }
}
static void BM_DeleteRowMultiplication(benchmark::State &state)
{
    int db_snps = state.range(0);
    serverInstance->GenData(1, db_snps);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        serverInstance->DeleteRowMultiplication(0);
    }
}

static void BM_SimilarityComputation(benchmark::State &state)
{
    int snps = state.range(0);
    int num_patients = state.range(1);

    vector<helib::Ctxt> p = vector<helib::Ctxt>();
    for (int i = 0; i < num_patients; i++)
    {
        p.push_back(serverInstance->Encrypt(0));
    }

    helib::Ctxt l1 = serverInstance->Encrypt(1);
    helib::Ctxt l2 = serverInstance->Encrypt(1);

    for (auto _ : state)
    {
        for (int j = 0; j < num_patients; j++)
        {
            for (int i = 0; i < snps; i++)
            {
                helib::Ctxt clone = l1;
                clone -= l2;
                clone.square();
                clone.cleanUp();
                p[j] += clone;
                benchmark::DoNotOptimize(p[j]);
                benchmark::DoNotOptimize(clone);
            }
        }
        benchmark::DoNotOptimize(p);
    }
}

BENCHMARK(BM_SimilarityComputation)->ArgsProduct({{100, 1000}, {1,2,3}})->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_GeneratePublicKeySwitch)->ArgsProduct({{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}})->Unit(benchmark::kSecond);
BENCHMARK(BM_SwithPublicKeySwitch)->ArgsProduct({{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}})->Unit(benchmark::kSecond);
BENCHMARK(BM_ParallelPRSQuery)->ArgsProduct({{1024, 4096, 16384}, benchmark::CreateRange(1, 16, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_ParallelMAFQuery)->ArgsProduct({{2, 4, 8, 16}, benchmark::CreateRange(1, 16, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_ParallelCountQuery)->ArgsProduct({{2, 4, 8, 16}, benchmark::CreateRange(1, 16, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_ParallelSimilarityQuery)->ArgsProduct({{1024, 4096, 16384}, {8}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_EncrpytCiphertext)->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_UpdateOneValue)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_UpdateOneRow)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_InsertRow)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_DeleteRowAddition)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_DeleteRowMultiplication)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK(BM_StorageCiphertext)->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK_MAIN();