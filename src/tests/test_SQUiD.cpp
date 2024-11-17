#include "../databases/plaintext_database.hpp"
#include "../databases/FHE_SIMD_database.hpp"
#include <gtest/gtest.h>

class SQUiDTest : public ::testing::Test
{
protected:
    // This is called before the first test
    static std::unique_ptr<FHESIMDDatabase> dbFHEInstance;
    static std::unique_ptr<PlaintextDatabase> dbInstance;

    static const uint32_t num_snp_cols = 3;
    static const uint32_t num_binary_pheno_cols = 1;
    static const uint32_t num_rows = 6;

    static void SetUpTestSuite()
    {
        dbFHEInstance = std::make_unique<FHESIMDDatabase>(constants::Large, false);
        dbInstance = std::make_unique<PlaintextDatabase>();

        std::cout << "Generating fake database..." << std::endl;
        const uint32_t seed = ::testing::UnitTest::GetInstance()->random_seed();

        dbFHEInstance->genData(num_rows, num_snp_cols, seed);
        dbInstance->genData(num_rows, num_snp_cols, seed);
    }

    // This is called after the last test
    static void TearDownTestSuite()
    {
    }

    // You can also define additional member variables or helper functions
    // that can be used in your tests.
};

std::unique_ptr<FHESIMDDatabase> SQUiDTest::dbFHEInstance = nullptr;
std::unique_ptr<PlaintextDatabase> SQUiDTest::dbInstance = nullptr; // Initialize to nullptr

TEST_F(SQUiDTest, PublicKeySwitch)
{
    Meta meta;
    meta(constants::Large);

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
    auto result_encrypted = SQUiDTest::dbFHEInstance->countQuery(1, query);
    auto result = SQUiDTest::dbFHEInstance->decrypt(result_encrypted)[0];

    auto true_count = SQUiDTest::dbInstance->countQuery(1, query);

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    cout << "Pred: " << result << endl;
    cout << "True: " << true_count << endl;

    ASSERT_EQ(true_count, result);

    const Meta &ref_Meta = SQUiDTest::dbFHEInstance->getMeta();

    helib::SecKey client_secret_key(ref_Meta.data->context);
    client_secret_key.GenSecKey();
    helib::PubKey client_public_key(client_secret_key);

    pair<vector<helib::DoubleCRT>, vector<helib::DoubleCRT>> ksk = client_public_key.genPublicKeySwitchingKey(ref_Meta.data->secretKey);

    helib::Ctxt ctxt(ref_Meta.data->publicKey);
    ctxt = result_encrypted;

    std::cout << "Noise before: " << ctxt.capacity() << std::endl;

    helib::Ctxt clone = ctxt;

    clone.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));

    std::cout << "Noise after: " << clone.capacity() << std::endl;

    helib::Ptxt<helib::BGV> new_plaintext_result(ref_Meta.data->context);
    client_secret_key.Decrypt(new_plaintext_result, clone);

    vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();

    long result2 = (long)poly_mod_result[0];
    ASSERT_EQ(result2, true_count);
}

TEST_F(SQUiDTest, ContextAndKeysSaveAndLoad)
{
    Meta meta;
    meta(constants::SmallNoSim);

    meta.data->saveToFile("test_context_and_keys.txt");
    Meta meta2;
    
    meta2("test_context_and_keys.txt", constants::SmallNoSim);

    ASSERT_EQ(meta.data->context, meta2.data->context);
}

// Add more tests as needed
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}