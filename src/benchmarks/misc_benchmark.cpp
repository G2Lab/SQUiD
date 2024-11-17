#include <iostream>
#include <vector>
#include <fstream>

#include <benchmark/benchmark.h>

#include "../databases/FHE_disk_database.hpp"
#include "../databases/plaintext_database.hpp"

static FHEDiskDatabase *dbFHEInstance;
static PlaintextDatabase *dbPlaintextInstance;
const uint32_t SEED = 0; 
const uint32_t NUM_ROWS = 49960;

static void DoSetup(const benchmark::State &state)
{
    static bool callSetup = true;
    if (callSetup)
    {
        dbFHEInstance = new FHEDiskDatabase(constants::Large, "", true);
        dbPlaintextInstance = new PlaintextDatabase();
        dbFHEInstance->printMeta();
        dbPlaintextInstance->genData(49960, 1024, SEED);
        dbPlaintextInstance->genBinaryPhenoData(1,SEED);

        dbFHEInstance->genData(49960, 16, SEED);

    }
    callSetup = false;
}

static void BM_EncrpytCiphertext(benchmark::State &state)
{
    std::vector<unsigned long> data = std::vector<unsigned long>(1000);

    for (auto _ : state)
    {
        auto result = dbFHEInstance->encrypt(data);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_StorageCiphertext(benchmark::State &state)
{
    std::vector<unsigned long> data = std::vector<unsigned long>(1);

    uint32_t result = 0;

    for (auto _ : state)
    {
        auto ctxt = dbFHEInstance->encrypt(data);
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
    Params param(constants::Large);

    uint32_t M = param.m;
    uint32_t P = param.p;
    uint32_t R = param.r;
    uint32_t BITS = param.qbits;
    
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                                .m(M)
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
        pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(owner_secret_key);
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
    state.counters["l"] = state.range(0);
}

static void BM_SwitchPublicKeySwitch(benchmark::State &state)
{
    Params param(constants::Large);
    uint32_t M = param.m;
    uint32_t P = param.p;
    uint32_t R = param.r;
    uint32_t BITS = param.qbits;
    
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                                .m(M)
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

    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(owner_secret_key);
    std::vector<helib::DoubleCRT>& firstVectorRef = ksk.first;
    std::vector<helib::DoubleCRT>& secondVectorRef = ksk.second;

    std::pair<std::vector<helib::DoubleCRT>&, std::vector<helib::DoubleCRT>&> kskRef(firstVectorRef, secondVectorRef);

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

        clone.PublicKeySwitch(kskRef);
        benchmark::DoNotOptimize(clone);
    }
    
    state.counters["l"] = state.range(0);
}

static void BM_ParallelSimilarityQuery(benchmark::State &state)
{
    vector<helib::Ctxt> d;
    for (int i = 0; i < state.range(0); i++)
    {
        d.push_back(dbFHEInstance->encrypt(0));
    }
    uint32_t targetSnp = 0;
    uint32_t threshold = 100;

    for (auto _ : state)
    {
        auto result = dbFHEInstance->similarityQueryP(targetSnp, d, threshold, state.range(1));

        state.PauseTiming();
        if (!result.first.isCorrect() || !result.second.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Threads"] = state.range(0);
    state.counters["Features"] = state.range(1);
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
        auto result = dbFHEInstance->PRSQueryP(query, state.range(1));

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Threads"] = state.range(1);
    state.counters["Features"] = state.range(0);
}

static void BM_ParallelMAFQuery(benchmark::State &state)
{
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }
    for (auto _ : state)
    {
        auto result = dbFHEInstance->MAFQueryP(0, query, state.range(1));

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }

        auto truth = dbPlaintextInstance->MAFQuery(0,1, query);
        auto result_decrypted = dbFHEInstance->decrypt(result);
        auto result_numerator = result_decrypted[0];
        auto result_denominator = result_decrypted[1];

        auto truth_numerator = truth % (2 * NUM_ROWS);
        auto truth_denominator = truth / (2 * NUM_ROWS);
        
        if (truth_numerator != result_numerator || truth_denominator != result_denominator)
        {
            std::cout << "Incorrect result " << std::endl;
            std::cout << "Numerator " << truth_numerator << " != " << result_numerator << std::endl;
            std::cout << "Denominator " << truth_denominator << " != " << result_denominator << std::endl;
        }

        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Threads"] = state.range(1);
    state.counters["Features"] = state.range(0);

}

static void BM_ParallelCountQuery(benchmark::State &state)
{
    vector<pair<uint32_t, uint32_t>> query = vector<pair<uint32_t, uint32_t>>();
    for (uint32_t i = 0; i < state.range(0); i++)
    {
        query.push_back(pair(i, 0));
    }
    for (auto _ : state)
    {
        auto result = dbFHEInstance->countQueryP(query, state.range(1));

        state.PauseTiming();
        if (!result.isCorrect())
        {
            std::cout << "ERROR EXCEEDED" << std::endl;
        }
        auto truth = dbPlaintextInstance->countQuery(1, query);
        auto result_decrypted = dbFHEInstance->decrypt(result)[0];
        if (truth != result_decrypted)
        {
            std::cout << "Incorrect result " << truth << " != " << result_decrypted << std::endl;
        }
        state.ResumeTiming();

        benchmark::DoNotOptimize(result);
    }
    state.counters["Threads"] = state.range(1);
    state.counters["Features"] = state.range(0);
}

static void BM_UpdateOneValue(benchmark::State &state)
{
    int db_snps = state.range(0);
    dbFHEInstance->genData(1, db_snps, SEED);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        dbFHEInstance->updateOneValue(0, 0, 0);
    }
}
static void BM_UpdateOneRow(benchmark::State &state)
{
    int db_snps = state.range(0);
    dbFHEInstance->genData(1, db_snps, SEED);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        dbFHEInstance->updateOneRow(0, vals);
    }
}

static void BM_InsertRow(benchmark::State &state)
{
    int db_snps = state.range(0);
    dbFHEInstance->genData(1, db_snps, SEED);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        dbFHEInstance->insertOneRow(vals);
    }
}
static void BM_DeleteRowAddition(benchmark::State &state)
{
    int db_snps = state.range(0);
    dbFHEInstance->genData(1, db_snps, SEED);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        dbFHEInstance->deleteRowAddition(0);
    }
}
static void BM_DeleteRowMultiplication(benchmark::State &state)
{
    int db_snps = state.range(0);
    dbFHEInstance->genData(1, db_snps, SEED);

    vector<uint32_t> vals = vector<uint32_t>(db_snps);
    for (int i = 0; i < db_snps; i++)
    {
        vals[i] = 0;
    }

    for (auto _ : state)
    {
        dbFHEInstance->deleteRowMultiplication(0);
    }
}

static void BM_SimilarityComputation(benchmark::State &state)
{
    int snps = state.range(0);
    int num_patients = state.range(1);

    vector<helib::Ctxt> p = vector<helib::Ctxt>();
    for (int i = 0; i < num_patients; i++)
    {
        p.push_back(dbFHEInstance->encrypt(0));
    }

    helib::Ctxt l1 = dbFHEInstance->encrypt(1);
    helib::Ctxt l2 = dbFHEInstance->encrypt(1);

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

BENCHMARK(BM_GeneratePublicKeySwitch)
   ->ArgsProduct({benchmark::CreateDenseRange(2, 20, 1)})
    ->Unit(benchmark::kSecond);

BENCHMARK(BM_SwitchPublicKeySwitch)
    ->ArgsProduct({benchmark::CreateDenseRange(2, 10, 1)})
    ->Unit(benchmark::kSecond);

BENCHMARK(BM_ParallelMAFQuery)->ArgsProduct({{2,4,8,16}, {1,2,4,8,16}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_ParallelCountQuery)->ArgsProduct({{2,4,8,16}, {1,2,4,8,16}})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_ParallelPRSQuery)->ArgsProduct({{1024, 4096, 16384}, benchmark::CreateRange(1, 16, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
BENCHMARK(BM_ParallelSimilarityQuery)->ArgsProduct({{1024, 4096, 16384}, {1,2,4,8,16}})->Unit(benchmark::kSecond)->Setup(DoSetup);

//BENCHMARK(BM_EncrpytCiphertext)->Unit(benchmark::kSecond)->Setup(DoSetup);
//BENCHMARK(BM_UpdateOneValue)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
//BENCHMARK(BM_UpdateOneRow)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
//BENCHMARK(BM_InsertRow)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
//BENCHMARK(BM_DeleteRowAddition)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
//BENCHMARK(BM_DeleteRowMultiplication)->ArgsProduct({benchmark::CreateRange(1, 1024, /*step=*/2)})->Unit(benchmark::kSecond)->Setup(DoSetup);
//BENCHMARK(BM_StorageCiphertext)->Unit(benchmark::kSecond)->Setup(DoSetup);

BENCHMARK_MAIN();
