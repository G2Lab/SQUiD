#include "aes_plaintext_database.hpp"
#include <stdexcept>

using namespace std;

void handleErrors(void)
{
    ERR_print_errors_fp(stderr);
    abort();
}

uint16_t decode_8bit(unsigned char *data)
{
    return (data[0]);
}
void *encode_8bit(uint16_t data, unsigned char *encoded)
{
    encoded[0] = data;
    return encoded;
}

int encrypt_128(unsigned char *plaintext, int plaintext_len, unsigned char *key,
                unsigned char *iv, unsigned char *ciphertext)
{
    EVP_CIPHER_CTX *ctx;

    int len;

    int ciphertext_len;

    /* Create and initialise the context */
    if (!(ctx = EVP_CIPHER_CTX_new()))
        handleErrors();

    /*
     * Initialise the encryption operation. IMPORTANT - ensure you use a key
     * and IV size appropriate for your cipher
     * In this example we are using 128-bit AES (i.e., a 128-bit key). The
     * IV size for *most* modes is the same as the block size. For AES this
     * is 128 bits
     */
    if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_128_cbc(), NULL, key, iv))
        handleErrors();

    /*
     * Provide the message to be encrypted, and obtain the encrypted output.
     * EVP_EncryptUpdate can be called multiple times if necessary
     */
    if (1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len))
        handleErrors();
    ciphertext_len = len;

    /*
     * Finalise the encryption. Further ciphertext bytes may be written at
     * this stage.
     */
    if (1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len))
        handleErrors();
    ciphertext_len += len;

    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);
    return ciphertext_len;
}

int decrypt_128(unsigned char *ciphertext, int ciphertext_len, unsigned char *key,
                unsigned char *iv, unsigned char *plaintext)
{
    EVP_CIPHER_CTX *ctx;

    int len;

    int plaintext_len;

    /* Create and initialise the context */
    if (!(ctx = EVP_CIPHER_CTX_new()))
        handleErrors();

    /*
     * Initialise the decryption operation. IMPORTANT - ensure you use a key
     * and IV size appropriate for your cipher
     * In this example, we are using 128-bit AES (i.e., a 128-bit key). The
     * IV size for *most* modes is the same as the block size. For AES, this
     * is 128 bits
     */
    if (1 != EVP_DecryptInit_ex(ctx, EVP_aes_128_cbc(), NULL, key, iv))
        handleErrors();

    /*
     * Provide the message to be decrypted, and obtain the plaintext output.
     * EVP_DecryptUpdate can be called multiple times if necessary.
     */
    if (1 != EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertext_len))
        handleErrors();
    plaintext_len = len;

    /*
     * Finalise the decryption. Further plaintext bytes may be written at
     * this stage.
     */
    if (1 != EVP_DecryptFinal_ex(ctx, plaintext + len, &len))
        handleErrors();
    plaintext_len += len;

    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);

    return plaintext_len;
}

PlaintextAESDatabase::PlaintextAESDatabase()
{
    initialise();
}

PlaintextAESDatabase::~PlaintextAESDatabase()
{
    delete[] key_128;
    delete[] iv;

    if (snp_data_set == true)
    {
        for (size_t i = 0; i < snp_data.size(); ++i)
        {
            for (size_t j = 0; j < snp_data[i].size(); ++j)
            {
                delete[] snp_data[i][j];
            }
        }
    }

    if (continuous_phenotype_data_set == true)
    {
        for (size_t i = 0; i < continuous_phenotype_data.size(); ++i)
        {
            for (size_t j = 0; j < continuous_phenotype_data[i].size(); ++j)
            {
                delete[] continuous_phenotype_data[i][j];
            }
        }
    }

    if (binary_phenotype_data_set == true)
    {
        for (size_t i = 0; i < binary_phenotype_data.size(); ++i)
        {
            for (size_t j = 0; j < binary_phenotype_data[i].size(); ++j)
            {
                delete[] binary_phenotype_data[i][j];
            }
        }
    }
}

void PlaintextAESDatabase::initialise()
{
    key_128 = new unsigned char[16]{0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
                                    0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35};

    /* A 128 bit IV */
    iv = new unsigned char[16]{0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
                               0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35};
}

void PlaintextAESDatabase::genData(uint32_t _num_rows, uint32_t _num_snp_cols, uint32_t seed)
{
    num_slots = 16;

    this->num_rows = _num_rows;
    this->num_snp_cols = _num_snp_cols;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 2);

    this->num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    snp_data.resize(_num_snp_cols, vector<unsigned char *>(num_compressed_rows, 0));

    for (size_t i = 0; i < snp_data.size(); ++i)
    {
        for (size_t j = 0; j < snp_data[i].size(); ++j)
        {
            snp_data[i][j] = new unsigned char[AES_BLOCK_SIZE];
            for (size_t k = 0; k < num_slots; ++k)
            {
                char val = dis(gen);
                encode_8bit(val, snp_data[i][j]);
            }
        }
    }

    snp_data_set = true;
}

void PlaintextAESDatabase::genContinuousPhenoData(uint32_t num_continuous_pheno_cols, uint32_t low, uint32_t high, uint32_t seed)
{
    this->num_continuous_pheno_cols = num_continuous_pheno_cols;

    continuous_phenotype_data.resize(num_continuous_pheno_cols, vector<unsigned char *>(num_compressed_rows, 0));

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(low, high);

    for (size_t i = 0; i < continuous_phenotype_data.size(); ++i)
    {
        for (size_t j = 0; j < continuous_phenotype_data[i].size(); ++j)
        {
            continuous_phenotype_data[i][j] = new unsigned char[AES_BLOCK_SIZE];
            for (size_t k = 0; k < num_slots; ++k)
            {
                char val = dis(gen);
                encode_8bit(val, continuous_phenotype_data[i][j]);
            }
        }
    }

    continuous_phenotype_data_set = true;
}

void PlaintextAESDatabase::genBinaryPhenoData(uint32_t num_binary_pheno_cols, uint32_t seed)
{
    this->num_binary_pheno_cols = num_binary_pheno_cols;

    binary_phenotype_data.resize(num_binary_pheno_cols, vector<unsigned char *>(num_rows, 0));

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 1);

    for (size_t i = 0; i < binary_phenotype_data.size(); ++i)
    {
        for (size_t j = 0; j < binary_phenotype_data[i].size(); ++j)
        {
            binary_phenotype_data[i][j] = new unsigned char[AES_BLOCK_SIZE];
            for (size_t k = 0; k < num_slots; ++k)
            {
                char val = dis(gen);
                encode_8bit(val, binary_phenotype_data[i][j]);
            }
        }
    }
}

void PlaintextAESDatabase::setData(vector<vector<int32_t>> &db)
{
    throw runtime_error("Not implemented");
}

void PlaintextAESDatabase::setBinaryPhenoData(vector<vector<int32_t>> &db)
{
    throw runtime_error("Not implemented");
}

void PlaintextAESDatabase::setData(string vcf_file)
{
    throw runtime_error("Not implemented");
}

unsigned char *PlaintextAESDatabase::getGenotype(uint32_t i, uint32_t j) const
{
    return snp_data[i][j];
}
unsigned char *PlaintextAESDatabase::getContinuousPheno(uint32_t i, uint32_t j) const
{
    return continuous_phenotype_data[i][j];
}
unsigned char *PlaintextAESDatabase::getBinaryPheno(uint32_t i, uint32_t j) const
{
    return binary_phenotype_data[i][j];
}

int32_t PlaintextAESDatabase::countQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const
{
    int count = 0;
    if (conjunctive)
    {
        for (size_t r = 0; r < num_compressed_rows; r++)
        {
            bool passes_filter[num_slots];
            memset(passes_filter, true, num_slots * sizeof(passes_filter));
            for (size_t c = 0; c < query.size(); c++)
            {
                int col = query[c].first;
                int value = query[c].second;
                unsigned char *decoded_result = getGenotype(col, r);
                for (size_t k = 0; k < num_slots; k++)
                {
                    int16_t spec_value = decode_8bit(decoded_result + k);
                    if (spec_value != value)
                    {
                        passes_filter[k] = false;
                    }
                }
            }
            for (size_t k = 0; k < num_slots; k++)
            {
                if (passes_filter[k])
                {
                    count += 1;
                }
            }
        }
    }
    else
    {
        for (size_t r = 0; r < num_compressed_rows; r++)
        {
            bool passes_filter[num_slots];
            memset(passes_filter, false, num_slots * sizeof(passes_filter));
            for (size_t c = 0; c < query.size(); c++)
            {
                int col = query[c].first;
                int value = query[c].second;
                unsigned char *decoded_result = getGenotype(col, r);
                for (size_t k = 0; k < num_slots; k++)
                {
                    int16_t spec_value = decode_8bit(decoded_result + k);
                    if (spec_value == value)
                    {
                        passes_filter[k] = true;
                    }
                }
            }
            for (size_t k = 0; k < num_slots; k++)
            {
                if (passes_filter[k])
                {
                    count += 1;
                }
            }
        }
    }
    return count;
}

int32_t PlaintextAESDatabase::MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const
{
    int alleles = 0;
    int patients = 0;
    if (conjunctive)
    {
        for (size_t r = 0; r < num_compressed_rows; r++)
        {
            bool passes_filter[num_slots];
            memset(passes_filter, true, num_slots * sizeof(passes_filter));
            for (size_t c = 0; c < query.size(); c++)
            {
                int col = query[c].first;
                int value = query[c].second;
                unsigned char *decoded_result = getGenotype(col, r);
                for (size_t k = 0; k < num_slots; k++)
                {
                    int16_t spec_value = decode_8bit(decoded_result + k);
                    if (static_cast<int>(decoded_result[k]) != value)
                    {
                        passes_filter[k] = false;
                    }
                }
            }
            unsigned char *decoded_target_snp = getGenotype(snp, r);

            for (size_t k = 0; k < num_slots; k++)
            {
                int16_t spec_value = decode_8bit(decoded_target_snp + k);
                if (passes_filter[k])
                {
                    alleles += spec_value;
                    patients += 2;
                }
            }
        }
    }
    else
    {
        for (size_t r = 0; r < num_compressed_rows; r++)
        {
            bool passes_filter[num_slots];
            memset(passes_filter, false, num_slots * sizeof(passes_filter));

            for (size_t c = 0; c < query.size(); c++)
            {
                int col = query[c].first;
                int value = query[c].second;
                unsigned char *decoded_result = getGenotype(col, r);
                for (size_t k = 0; k < num_slots; k++)
                {
                    int16_t spec_value = decode_8bit(decoded_result + k);
                    if (spec_value == value)
                    {
                        passes_filter[k] = true;
                    }
                }
            }
            unsigned char *decoded_target_snp = getGenotype(snp, r);

            for (size_t k = 0; k < num_slots; k++)
            {
                int16_t spec_value = decode_8bit(decoded_target_snp + k);
                if (passes_filter[k])
                {
                    alleles += spec_value;
                    patients += 2;
                }
            }
        }
    }
    return patients + alleles * (2 * num_rows);
}

vector<int32_t> PlaintextAESDatabase::PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params) const
{
    vector<int32_t> scores = vector<int32_t>(num_rows);

    for (size_t r = 0; r < num_compressed_rows; r++)
    {
        int score[num_slots];
        memset(score, 0, num_slots * sizeof(int));

        for (size_t c = 0; c < prs_params.size(); c++)
        {
            int col = prs_params[c].first;
            int value = prs_params[c].second;
            unsigned char *decoded_result = getGenotype(col, r);
            for (size_t k = 0; k < num_slots; k++)
            {
                uint16_t val = decode_8bit(decoded_result + k);
                score[k] += value * val;
            }
        }
        for (size_t k = 0; k < num_slots; k++)
        {
            if ((num_slots * r + k) < num_rows)
            {
                scores.at(num_slots * r + k) = score[k];
            }
        }
    }
    return scores;
}

pair<int32_t, int32_t> PlaintextAESDatabase::similarityQuery(uint32_t target_column, vector<unsigned char *> &d, uint32_t threshold) const
{
    int32_t yes = 0;
    int32_t no = 0;

    for (size_t r = 0; r < num_compressed_rows; r++)
    {
        int score[num_slots];
        memset(score, 0, num_slots * sizeof(int));

        for (size_t c = 0; c < d.size(); c++)
        {
            int16_t d_val = decode_8bit(d[c]);
            unsigned char *decoded_result = getGenotype(c, r);
            for (size_t k = 0; k < num_slots; k++)
            {
                int16_t val = decode_8bit(decoded_result + k);
                score[k] += pow(val - d_val, 2);
            }
        }
        unsigned char *decoded_target_column = getBinaryPheno(target_column, r);
        for (size_t k = 0; k < num_slots; k++)
        {
            int16_t val = decode_8bit(decoded_target_column + k);
            if (score[k] < threshold)
            {
                yes += val;
                no += 1 - val;
            }
        }
    }
    return pair(yes, no);
}

int32_t PlaintextAESDatabase::countingRangeQuery(uint32_t lower, uint32_t upper, uint32_t phenotype_index)
{
    throw runtime_error("Not implemented");
}

pair<int32_t, int32_t> PlaintextAESDatabase::MAFRangeQuery(uint32_t snp, uint32_t lower, uint32_t upper, uint32_t phenotype_index)
{
    throw runtime_error("Not implemented");
}

void PlaintextAESDatabase::printData(bool with_headers) const
{
    if (with_headers)
    {
        for (const auto &header : column_headers)
        {
            cout << header << " ";
        }
        cout << endl;
    }

    for (size_t i = 0; i < snp_data[0].size(); ++i)
    {
        for (size_t j = 0; j < snp_data.size(); ++j)
        {
            cout << snp_data[j][i] << " ";
        }
        cout << endl;
    }
}

void PlaintextAESDatabase::printMeta() const
{
    cout << "Number of rows: " << num_rows << endl;
    cout << "Number of snp columns: " << num_snp_cols << endl;
    cout << "Number of binary pheno columns: " << num_binary_pheno_cols << endl;
    cout << "Number of continuous pheno columns: " << num_continuous_pheno_cols << endl;
    cout << "Number of compressed rows: " << num_compressed_rows << endl;
}

void PlaintextAESDatabase::setGenotype(unsigned char *data, uint32_t column, uint32_t compressed_row_index)
{
    for (size_t k = 0; k < num_slots; k++)
    {
        for (size_t i = 0; i < 2; i++)
        {
            snp_data[column][compressed_row_index][k * 2 + i] = data[k * 2 + i];
        }
    }
}