#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>

#include "server.hpp"
#include "globals.hpp"
#include "tools.hpp"

class SQUiDTest : public ::testing::Test
{
protected:
    // This is called before the first test
    static std::unique_ptr<Server> serverInstance;
    static std::vector<std::vector<uint32_t>> *fake_db; // Change to pointer

    static const int num_cols = 3;
    static const int num_rows = 100;

    static void SetUpTestSuite()
    {
        serverInstance = std::make_unique<Server>(constants::P131, true);

        std::cout << "Generating fake database..." << std::endl;

        // Allocate memory for fake_db
        fake_db = new std::vector<std::vector<uint32_t>>(num_cols, std::vector<uint32_t>(num_rows, 0));

        // Initialize fake_db with random values
        for (int i = 0; i < num_cols; i++)
        {
            for (int j = 0; j < num_rows; j++)
            {
                (*fake_db)[i][j] = rand() % 2;
            }
        }

        serverInstance->SetData(*fake_db);
    }

    // This is called after the last test
    static void TearDownTestSuite()
    {
        // Deallocate memory for fake_db
        delete fake_db;
    }

    // You can also define additional member variables or helper functions
    // that can be used in your tests.
};

std::unique_ptr<Server> SQUiDTest::serverInstance = nullptr;
std::vector<std::vector<uint32_t>> *SQUiDTest::fake_db = nullptr; // Initialize to nullptr

TEST_F(SQUiDTest, CountingQueryAnd)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = SQUiDTest::serverInstance->CountQuery(1, query);
    auto result = SQUiDTest::serverInstance->Decrypt(result_encrypted)[0];

    int true_count = 0;
    for (int i = 0; i < num_rows; i++)
    {
        if ((*fake_db)[0][i] == 0 && (*fake_db)[1][i] == 1)
        {
            true_count++;
        }
    }

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(SQUiDTest, CountingQueryOr)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = SQUiDTest::serverInstance->CountQuery(0, query);
    auto result = SQUiDTest::serverInstance->Decrypt(result_encrypted)[0];

    int true_count = 0;
    for (int i = 0; i < num_rows; i++)
    {
        if ((*fake_db)[0][i] == 0 || (*fake_db)[1][i] == 1)
        {
            true_count++;
        }
    }

    cout << "Running Counting query (snp 0 = 0 or snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);
}

TEST_F(SQUiDTest, MAFQuery)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 1)};
    auto result_encrypted = SQUiDTest::serverInstance->MAFQuery(2, 1, query);
    auto result = SQUiDTest::serverInstance->Decrypt(result_encrypted);
    auto nom = result[0];
    auto dom = result[1];
    auto MAF = (double)min(nom, dom - nom) / dom;

    int true_count = 0;
    int passing_rows = 0;
    for (int i = 0; i < num_rows; i++)
    {
        if ((*fake_db)[0][i] == 1)
        {
            true_count += (*fake_db)[2][i];
            passing_rows++;
        }
    }

    passing_rows *= 2;

    double true_maf = (double)min(true_count, passing_rows - true_count) / passing_rows;

    cout << "Running MAF query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << MAF << endl;
    cout << "True: " << true_maf << endl;

    ASSERT_EQ(true_maf, MAF);
}

TEST_F(SQUiDTest, PRSQuery)
{
    vector<pair<uint32_t, int>> query;
    query = vector<pair<uint32_t, int>>{pair(0, 2), pair(1, 3), pair(2, 9)};
    auto result_encrypted = SQUiDTest::serverInstance->PRSQuery(query);
    auto result = SQUiDTest::serverInstance->Decrypt(result_encrypted[0]);

    std::vector<int> true_result(SQUiDTest::num_rows, 0);
    for (int i = 0; i < SQUiDTest::num_rows; i++)
    {
        true_result[i] = 2 * (*fake_db)[0][i] + 3 * (*fake_db)[1][i] + 9 * (*fake_db)[2][i];
    }

    for (int i = 0; i < 100; i++)
    {
        ASSERT_EQ(true_result[i], result[i]);
    }
}

TEST_F(SQUiDTest, SimilarityQuery)
{
    vector<helib::Ctxt> d = vector<helib::Ctxt>();

    for (int i = 0; i < 2; i++)
    {
        d.push_back(SQUiDTest::serverInstance->Encrypt(2));
    }
    int threshold = 2;

    auto result_encrypted = SQUiDTest::serverInstance->SimilarityQuery(2, d, threshold);
    auto with = SQUiDTest::serverInstance->Decrypt(result_encrypted.first)[0];
    auto without = SQUiDTest::serverInstance->Decrypt(result_encrypted.second)[0];

    int true_with = 0;
    int true_without = 0;
    for (int i = 0; i < num_rows; i++)
    {

        if (pow((*fake_db)[0][i] - 2, 2) + pow((*fake_db)[1][i] - 2, 2) <= threshold)
        {
            if ((*fake_db)[2][i] == 2)
            {
                true_with++;
            }
            else
            {
                true_without++;
            }
        }
    }

    cout << "Running similarity query (d: snp 0 = 2 and snp 1 = 2, target = 2, threshold = 8)" << endl;
    cout << "Count with target:   " << with << endl;
    cout << "Count without target:" << without << endl;
    cout << "True with: " << true_with << endl;
    cout << "True without: " << true_without << endl;

    ASSERT_EQ(true_with, with);
}

TEST_F(SQUiDTest, PublicKeySwitch)
{
    Meta meta;
    meta(constants::P131);

    helib::SecKey owner_secret_key(meta.data->context);
    owner_secret_key.GenSecKey();
    helib::PubKey owner_public_key(owner_secret_key);

    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);

    auto ksk = client_public_key.genPublicKeySwitchingKey(owner_secret_key);

    helib::Ptxt<helib::BGV> ptxt(meta.data->context);
    int num_slots = 10;
    vector<long> original_values = vector<long>(num_slots);
    for (int i = 0; i < num_slots; i++)
    {
        ptxt[i] = i;
        original_values[i] = i;
    }

    helib::Ctxt ctxt(owner_public_key);
    owner_public_key.Encrypt(ctxt, ptxt);

    std::cout << "Noise before: " << ctxt.capacity() << std::endl;

    helib::Ctxt clone = ctxt;

    clone.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));

    std::cout << "Noise after: " << clone.capacity() << std::endl;

    helib::Ptxt<helib::BGV> new_plaintext_result(meta.data->context);
    client_secret_key.Decrypt(new_plaintext_result, clone);

    vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();

    vector<long> result = vector<long>(num_slots);

    for (uint32_t i = 0; i < num_slots; i++)
    {
        result[i] = (long)poly_mod_result[i];
    }
    for (int i = 0; i < num_slots; i++)
    {
        ASSERT_EQ(result[i], original_values[i]);
    }
}
TEST_F(SQUiDTest, CountAndKeySwitch)
{
    vector<pair<uint32_t, uint32_t>> query;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    auto result_encrypted = SQUiDTest::serverInstance->CountQuery(1, query);
    auto result = SQUiDTest::serverInstance->Decrypt(result_encrypted)[0];

    int true_count = 0;
    for (int i = 0; i < num_rows; i++)
    {
        if ((*fake_db)[0][i] == 0 && (*fake_db)[1][i] == 1)
        {
            true_count++;
        }
    }

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);

    Meta meta;
    meta(constants::P131);

    helib::SecKey client_secret_key(meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);

    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(SQUiDTest::serverInstance->GetMeta().data->secretKey);

    helib::Ctxt ctxt(SQUiDTest::serverInstance->GetMeta().data->publicKey);
    ctxt = result_encrypted;

    std::cout << "Noise before: " << ctxt.capacity() << std::endl;

    helib::Ctxt clone = ctxt;

    clone.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));

    std::cout << "Noise after: " << clone.capacity() << std::endl;

    helib::Ptxt<helib::BGV> new_plaintext_result(meta.data->context);
    client_secret_key.Decrypt(new_plaintext_result, clone);

    vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();

    vector<long> result2 = vector<long>(1);

    for (uint32_t i = 0; i < 100; i++)
    {
        result2[i] = (long)poly_mod_result[i];
    }
    for (int i = 0; i < 100; i++)
    {
        ASSERT_EQ(result2[i], true_count);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}