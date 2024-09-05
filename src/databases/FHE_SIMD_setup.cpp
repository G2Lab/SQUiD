#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

void FHESIMDDatabase::genData(uint32_t num_rows, uint32_t num_snp_cols, uint32_t seed)
{
    this->num_rows = num_rows;
    this->num_snp_cols = num_snp_cols;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 2);

    this->num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;
    
    this->snp_data = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_snp_cols; i++)
    {
        uint32_t rows_set = 0;
        vector<helib::Ctxt> snp_data_column = vector<helib::Ctxt>();
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = dis(gen);
                rows_set += 1;
            }

            helib::Ctxt ctxt = encrypt(ptxt);

            snp_data_column.push_back(ctxt);
        }
        this->snp_data.push_back(snp_data_column);
    }
    this->snp_data_set = true;
}

void FHESIMDDatabase::genContinuousPhenoData(uint32_t num_continuous_pheno_cols, uint32_t low, uint32_t high, uint32_t seed)
{
    this->num_continuous_pheno_cols = num_continuous_pheno_cols;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(low, high);

    continuous_phenotype_data = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_continuous_pheno_cols; i++){
        vector<helib::Ctxt> phenotype_data_column = vector<helib::Ctxt>();
        uint32_t rows_set = 0;
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ptxt<helib::BGV> ptxt(meta.data->context);
            helib::Ctxt ctxt(meta.data->publicKey);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = dis(gen);
                rows_set += 1;
            }

            meta.data->publicKey.Encrypt(ctxt, ptxt.getPolyRepr());

            phenotype_data_column.push_back(ctxt);
        }
        continuous_phenotype_data.push_back(phenotype_data_column);
    }

    continuous_phenotype_data_set = true;
}

void FHESIMDDatabase::genBinaryPhenoData(uint32_t num_binary_pheno_cols, uint32_t seed)
{
    this->num_binary_pheno_cols = num_binary_pheno_cols;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 1);

    this->binary_phenotype_data = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_binary_pheno_cols; i++){
        vector<helib::Ctxt> phenotype_data_column = vector<helib::Ctxt>();
        uint32_t rows_set = 0;
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ptxt<helib::BGV> ptxt(meta.data->context);
            helib::Ctxt ctxt(meta.data->publicKey);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = dis(gen);
                rows_set += 1;
            }

            meta.data->publicKey.Encrypt(ctxt, ptxt.getPolyRepr());

            phenotype_data_column.push_back(ctxt);
        }
        this->binary_phenotype_data.push_back(phenotype_data_column);
    }
    this->binary_phenotype_data_set = true;
}

void FHESIMDDatabase::setData(vector<vector<int32_t>> &db)
{
    num_snp_cols = db.size();
    if (num_snp_cols == 0)
    {
        throw invalid_argument("ERROR: DB has zero columns! THIS DOES NOT WORK!");
    }

    num_rows = db[0].size();

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    snp_data = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_snp_cols; i++)
    {
        vector<helib::Ctxt> cipher_vector = vector<helib::Ctxt>();
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {

            vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = db[i][j * num_slots + k];
            }

            helib::Ctxt ctxt = encrypt(ptxt);

            cipher_vector.push_back(ctxt);
        }
        snp_data.push_back(cipher_vector);
    }

    snp_data_set = true;
}

void FHESIMDDatabase::setBinaryPhenoData(vector<vector<int32_t>> &db)
{
    num_binary_pheno_cols = db.size();
    if (num_binary_pheno_cols == 0)
    {
        throw invalid_argument("ERROR: DB has zero columns! THIS DOES NOT WORK!");
    }

    binary_phenotype_data = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_binary_pheno_cols; i++)
    {
        vector<helib::Ctxt> pheno_column = vector<helib::Ctxt>();
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {

            vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = db[i][j * num_slots + k];
            }

            helib::Ctxt ctxt = encrypt(ptxt);

            pheno_column.push_back(ctxt);
        }
        binary_phenotype_data.push_back(pheno_column);
    }

    binary_phenotype_data_set = true;
}

void FHESIMDDatabase::setData(string vcf_file)
{
    std::ifstream file(vcf_file);
    if (!file.is_open())
    {
        std::cout << "Error opening file: " << vcf_file << std::endl;
        return;
    }

    uint32_t col_counter = 0;
    uint32_t row_counter = 0;
    uint32_t delimiter_counter = 0;

    std::vector<std::vector<int32_t>> matrix;
    std::string line;
    while (std::getline(file, line))
    {
        // Skip header lines starting with "#"
        if (line[0] == '#')
            continue;

        std::vector<int32_t> row;

        col_counter += 1;

        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, '\t'))
        {
            if (delimiter_counter == 2)
            { // When we are reading in the snp name
                column_headers.push_back(token);
            }

            delimiter_counter += 1;
            if (token.find(".") != std::string::npos)
            {
                row.push_back(-1);
                row_counter += 1;
                continue;
            }
            if (token == "1/1" || token == "1|1")
            {
                row.push_back(2);
                row_counter += 1;
            }
            if (token == "0/1" || token == "0|1")
            {
                row.push_back(1);
                row_counter += 1;
            }
            if (token == "1|0")
            {
                row.push_back(1);
                row_counter += 1;
            }
            if (token == "0/0" || token == "0|0")
            {
                row.push_back(0);
                row_counter += 1;
            }
        }
        num_rows = row_counter;
        row_counter = 0;
        delimiter_counter = 0;

        matrix.push_back(row);
    }
    num_snp_cols = col_counter;

    setData(matrix);
    file.close();

    snp_data_set = true;
}

void FHESIMDDatabase::setColumnHeaders(vector<string> &headers)
{
    column_headers = vector<string>();
    for (uint32_t i = 0; i < headers.size(); i++)
    {
        column_headers.push_back(headers[i]);
    }
}