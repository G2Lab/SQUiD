#include "../databases/plaintext_database.hpp"
#include "../databases/FHE_SIMD_database.hpp"
#include <gtest/gtest.h>

long aggregate(const std::vector<long> &v, long entries)
{
    long sum = 0;
    for (long i = 0; i < entries; i++)
    {
        sum += v[i];
    }
    return sum;
}


class FHESIMDDatabaseTest : public ::testing::Test
{
protected:
    // This is called before the first test
    static std::unique_ptr<FHESIMDDatabase> dbFHEInstance;
    static std::unique_ptr<PlaintextDatabase> dbInstance;

    static const uint32_t num_snp_cols = 3;
    static const uint32_t num_binary_pheno_cols = 1;
    static const uint32_t num_rows = 600;
    static const uint32_t num_threads = 2;
    static const uint32_t low = 0;
    static const uint32_t high = 10;

    static void SetUpTestSuite()
    {
        dbFHEInstance = std::make_unique<FHESIMDDatabase>(constants::Large, true);
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



std::unique_ptr<FHESIMDDatabase> FHESIMDDatabaseTest::dbFHEInstance = nullptr;
std::unique_ptr<PlaintextDatabase> FHESIMDDatabaseTest::dbInstance = nullptr; // Initialize to nullptr

TEST_F(FHESIMDDatabaseTest, SimilarityQuery)
{
    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    vector<int32_t> d_plain = vector<int32_t>();

    for (int i = 0; i < 2; i++)
    {
        d.push_back(FHESIMDDatabaseTest::dbFHEInstance->encrypt(2));
        d_plain.push_back(2);
    }

    int threshold = 3;
    int target_column = 0;

    auto result_encrypted = FHESIMDDatabaseTest::dbFHEInstance->similarityQuery(target_column, d, threshold);

    auto with = FHESIMDDatabaseTest::dbFHEInstance->decrypt(result_encrypted.first);
    auto without = FHESIMDDatabaseTest::dbFHEInstance->decrypt(result_encrypted.second);

    auto with_long = aggregate(with, result_encrypted.first.nAggregates);
    auto without_long = aggregate(without, result_encrypted.second.nAggregates);

    auto true_result = FHESIMDDatabaseTest::dbInstance->similarityQuery(target_column, d_plain, threshold);

    int true_with = true_result.first;
    int true_without = true_result.second;

    cout << "Running similarity query (d: snp 0 = 2 and snp 1 = 2, target = 0, threshold = 2)" << endl;
    cout << "Count with target:   " << with_long << endl;
    cout << "Count without target:" << without_long << endl;
    cout << "True with: " << true_with << endl;
    cout << "True without: " << true_without << endl;

    ASSERT_EQ(true_with, with_long);
    ASSERT_EQ(true_without, without_long);
}

TEST_F(FHESIMDDatabaseTest, SimilarityQueryParallel)
{
    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    vector<int32_t> d_plain = vector<int32_t>();

    for (int i = 0; i < 2; i++)
    {
        d.push_back(FHESIMDDatabaseTest::dbFHEInstance->encrypt(2));
        d_plain.push_back(2);
    }

    int threshold = 3;
    int target_column = 0;

    auto result_encrypted = FHESIMDDatabaseTest::dbFHEInstance->similarityQueryP(target_column, d, threshold, num_threads);

    std::cout << "Capacity: " << result_encrypted.first.capacity() << std::endl;
    std::cout << "Capacity: " << result_encrypted.second.capacity() << std::endl;

    auto with = FHESIMDDatabaseTest::dbFHEInstance->decrypt(result_encrypted.first);
    auto without = FHESIMDDatabaseTest::dbFHEInstance->decrypt(result_encrypted.second);

    auto true_result = FHESIMDDatabaseTest::dbInstance->similarityQuery(target_column, d_plain, threshold);

    auto with_long = aggregate(with, result_encrypted.first.nAggregates);
    auto without_long = aggregate(without, result_encrypted.second.nAggregates);

    int true_with = true_result.first;
    int true_without = true_result.second;

    cout << "Running similarity query (d: snp 0 = 2 and snp 1 = 2, target = 0, threshold = 2)" << endl;
    cout << "Count with target:   " << with_long << endl;
    cout << "Count without target:" << without_long << endl;
    cout << "True with: " << true_with << endl;
    cout << "True without: " << true_without << endl;

    ASSERT_EQ(true_with, with_long);
    ASSERT_EQ(true_without, without_long);
}

TEST_F(FHESIMDDatabaseTest, CountRangeQuery)
{
    uint32_t low_query = 2;
    uint32_t high_query = 5;
    auto result_encrypted = FHESIMDDatabaseTest::dbFHEInstance->countingRangeQuery(low_query, high_query, 0);
    auto result = FHESIMDDatabaseTest::dbFHEInstance->decrypt(result_encrypted);

    auto result_long = aggregate(result, result_encrypted.nAggregates);

    uint32_t true_count = FHESIMDDatabaseTest::dbInstance->countingRangeQuery(low_query, high_query, 0);

    cout << "Running Counting query (pheno 0 in [2, 5])" << endl;
    cout << "Pred: " << result_long << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result_long);
}

TEST_F(FHESIMDDatabaseTest, MAFRangeQuery)
{
    uint32_t low_query = 2;
    uint32_t high_query = 5;
    auto result_encrypted = FHESIMDDatabaseTest::dbFHEInstance->MAFRangeQuery(2, low_query, high_query, 0);
    auto result_num = FHESIMDDatabaseTest::dbFHEInstance->decrypt(result_encrypted.first);
    auto result_dom = FHESIMDDatabaseTest::dbFHEInstance->decrypt(result_encrypted.second);
    auto true_result = FHESIMDDatabaseTest::dbInstance->MAFRangeQuery(2, low_query, high_query, 0);

    auto num = aggregate(result_num, result_encrypted.first.nAggregates);
    auto dom = aggregate(result_dom, result_encrypted.second.nAggregates);

    uint32_t true_num = true_result.second;
    uint32_t true_dom = true_result.first;


    cout << "Running MAF query (pheno 0 in [2, 5], for snp 2)" << endl;
    cout << "Pred: " << num << "/" << dom << endl;
    cout << "True: " << true_num << "/" << true_dom << endl;

    ASSERT_EQ(true_num, num);
    ASSERT_EQ(true_dom, dom);
}

// Add more tests as needed
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}