#include "../databases/plaintext_database.hpp"
#include "../databases/FHE_SIMD_database.hpp"
#include <gtest/gtest.h>

class FHESIMDDatabaseTestNoComp : public ::testing::Test
{
protected:
    // This is called before the first test
    static std::unique_ptr<FHESIMDDatabase> dbFHEInstance;
    static std::unique_ptr<PlaintextDatabase> dbInstance;

    static const uint32_t num_snp_cols = 3;
    static const uint32_t num_binary_pheno_cols = 1;
    static const uint32_t num_rows = 6;
    static const uint32_t num_threads = 2;
    static const uint32_t low = 0;
    static const uint32_t high = 10;

    static void SetUpTestSuite()
    {
        dbFHEInstance = std::make_unique<FHESIMDDatabase>(constants::Test, false);
        dbInstance = std::make_unique<PlaintextDatabase>();

        std::cout << "Generating fake database..." << std::endl;
        const uint32_t seed = ::testing::UnitTest::GetInstance()->random_seed();

        dbFHEInstance->genData(num_rows, num_snp_cols, seed);
        dbInstance->genData(num_rows, num_snp_cols, seed);

        dbFHEInstance->genBinaryPhenoData(num_binary_pheno_cols, seed);
        dbInstance->genBinaryPhenoData(num_binary_pheno_cols, seed);

        dbFHEInstance->genContinuousPhenoData(1, low, high, seed);
        dbInstance->genContinuousPhenoData(1, low, high, seed);
    }

    // This is called after the last test
    static void TearDownTestSuite()
    {
    }

    // You can also define additional member variables or helper functions
    // that can be used in your tests.
};

std::unique_ptr<FHESIMDDatabase> FHESIMDDatabaseTestNoComp::dbFHEInstance = nullptr;
std::unique_ptr<PlaintextDatabase> FHESIMDDatabaseTestNoComp::dbInstance = nullptr; // Initialize to nullptr

TEST_F(FHESIMDDatabaseTestNoComp, CountingQueryAnd)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->countQuery(1, query);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHESIMDDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHESIMDDatabaseTestNoComp, CountingQueryOr)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->countQuery(1, query);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHESIMDDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 or snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHESIMDDatabaseTestNoComp, CountingQueryOr2)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 2)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->countQuery(1, query);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHESIMDDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 or snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHESIMDDatabaseTestNoComp, MAFQuery)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 1)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->MAFQuery(2, 1, query);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted);
    auto true_result = FHESIMDDatabaseTestNoComp::dbInstance->MAFQuery(2, 1, query);

    auto num = result[0];
    auto dom = result[1];

    uint32_t true_num = true_result / (2 * num_rows);
    uint32_t true_dom = true_result % (2 * num_rows);


    cout << "Running MAF query (snp 0 = 1)" << endl;
    cout << "Pred: " << num << "/" << dom << endl;
    cout << "True: " << true_num << "/" << true_dom << endl;

    ASSERT_EQ(true_num, num);
    ASSERT_EQ(true_dom, dom);
}

TEST_F(FHESIMDDatabaseTestNoComp, PRSQuery)
{
    vector<pair<uint32_t, int>> query;
    query = vector<pair<uint32_t, int>>{pair(0, 2), pair(1, 3), pair(2, 1)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->PRSQuery(query);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted[0]);

    auto true_result = FHESIMDDatabaseTestNoComp::dbInstance->PRSQuery(query);

    for (uint32_t i = 0; i < num_rows; i++)
    {
        ASSERT_EQ(true_result[i], result[i]);
    }
}

TEST_F(FHESIMDDatabaseTestNoComp, CountingQueryParallel)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->countQueryP(query, num_threads);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHESIMDDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHESIMDDatabaseTestNoComp, MAFQueryParallel)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->MAFQueryP(0, query, num_threads);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted);
    auto true_result = FHESIMDDatabaseTestNoComp::dbInstance->MAFQuery(0, 1, query);

    auto num = result[0];
    auto dom = result[1];

    uint32_t true_num = true_result / (2 * num_rows);
    uint32_t true_dom = true_result % (2 * num_rows);


    cout << "Running MAF query (snp 0 = 1)" << endl;
    cout << "Pred: " << num << "/" << dom << endl;
    cout << "True: " << true_num << "/" << true_dom << endl;

    ASSERT_EQ(true_num, num);
    ASSERT_EQ(true_dom, dom);
}

TEST_F(FHESIMDDatabaseTestNoComp, PRSQueryParallel)
{
    vector<pair<uint32_t, int>> query;
    query = vector<pair<uint32_t, int>>{pair(0, 2), pair(1, 3), pair(2, 1)};
    auto result_encrypted = FHESIMDDatabaseTestNoComp::dbFHEInstance->PRSQueryP(query, num_threads);
    auto result = FHESIMDDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted);

    auto true_result = FHESIMDDatabaseTestNoComp::dbInstance->PRSQuery(query);

    for (uint32_t i = 0; i < num_rows; i++)
    {
        ASSERT_EQ(true_result[i], result[i]);
    }
}


// Add more tests as needed
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}