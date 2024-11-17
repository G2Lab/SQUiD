#include "../databases/plaintext_database.hpp"
#include "../databases/FHE_disk_database.hpp"
#include <gtest/gtest.h>

class FHEDIskDatabaseTestNoComp : public ::testing::Test
{
protected:
    // This is called before the first test
    static std::unique_ptr<FHEDiskDatabase> dbFHEInstance;
    static std::unique_ptr<PlaintextDatabase> dbInstance;

    static const uint32_t num_snp_cols = 3;
    static const uint32_t num_binary_pheno_cols = 1;
    static const uint32_t num_rows = 5;
    static const uint32_t num_threads = 2;
    static const uint32_t low = 0;
    static const uint32_t high = 10;

    static void SetUpTestSuite()
    {
        dbFHEInstance = std::make_unique<FHEDiskDatabase>(constants::Test, false);
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

std::unique_ptr<FHEDiskDatabase> FHEDIskDatabaseTestNoComp::dbFHEInstance = nullptr;
std::unique_ptr<PlaintextDatabase> FHEDIskDatabaseTestNoComp::dbInstance = nullptr; // Initialize to nullptr

TEST_F(FHEDIskDatabaseTestNoComp, CountingQueryAnd)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->countQuery(1, query);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHEDIskDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHEDIskDatabaseTestNoComp, CountingQueryOr)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->countQuery(1, query);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHEDIskDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 or snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHEDIskDatabaseTestNoComp, CountingQueryOr2)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 2)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->countQuery(1, query);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHEDIskDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 or snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHEDIskDatabaseTestNoComp, PRSQuery)
{
    vector<pair<uint32_t, int>> query;
    query = vector<pair<uint32_t, int>>{pair(0, 2), pair(1, 3), pair(2, 1)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->PRSQuery(query);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted[0]);

    auto true_result = FHEDIskDatabaseTestNoComp::dbInstance->PRSQuery(query);

    for (uint32_t i = 0; i < num_rows; i++)
    {
        ASSERT_EQ(true_result[i], result[i]);
    }
}

TEST_F(FHEDIskDatabaseTestNoComp, CountingQueryParallel)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->countQueryP(query, num_threads);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted)[0];

    uint32_t true_count = FHEDIskDatabaseTestNoComp::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(FHEDIskDatabaseTestNoComp, MAFQuery)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 1)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->MAFQuery(2, 1, query);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted);
    auto true_result = FHEDIskDatabaseTestNoComp::dbInstance->MAFQuery(2, 1, query);

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

TEST_F(FHEDIskDatabaseTestNoComp, MAFQueryParallel)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 1)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->MAFQueryP(2, query, num_threads);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted);
    auto true_result = FHEDIskDatabaseTestNoComp::dbInstance->MAFQuery(2, 1, query);

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

TEST_F(FHEDIskDatabaseTestNoComp, PRSQueryParallel)
{
    vector<pair<uint32_t, int>> query;
    query = vector<pair<uint32_t, int>>{pair(0, 2), pair(1, 3), pair(2, 1)};
    auto result_encrypted = FHEDIskDatabaseTestNoComp::dbFHEInstance->PRSQueryP(query, num_threads);
    auto result = FHEDIskDatabaseTestNoComp::dbFHEInstance->decrypt(result_encrypted);

    auto true_result = FHEDIskDatabaseTestNoComp::dbInstance->PRSQuery(query);

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