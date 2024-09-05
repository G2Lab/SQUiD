#include "../databases/plaintext_database.hpp"
#include <gtest/gtest.h>

TEST(PlaintextDatabaseTest, GenerateData) {
    PlaintextDatabase db;
    db.genData(10, 5, 0);

    ASSERT_EQ(db.num_rows, 10);
    ASSERT_EQ(db.num_snp_cols, 5);
}

TEST(PlaintextDatabaseTest, SimilarityQuery) {
    PlaintextDatabase db;
    vector<vector<int32_t>> fake_db = vector<vector<int32_t>>{
        {0, 0, 0, 1, 1, 1, 2, 2, 2, 0},
        {0, 1, 2, 0, 1, 2, 0, 1, 2, 1}};
    
    vector<vector<int32_t>> fake_binary_pheno_db = vector<vector<int32_t>>{
        {0, 0, 0, 0, 1, 1, 0, 1, 0, 0}};

    db.setData(fake_db);
    db.setBinaryPhenoData(fake_binary_pheno_db);

    vector<int32_t> d = vector<int32_t>();
    
    d.push_back(2);
    d.push_back(2);
    
    auto result_sim = db.similarityQuery(0, d, 4);
    ASSERT_EQ(result_sim.first, 3);
    ASSERT_EQ(result_sim.second, 1);
}

// Add more tests as needed
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
