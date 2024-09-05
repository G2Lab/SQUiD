#include "FHE_disk_database.hpp"

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------

//                                                         HELPER FUNCTIONS

// ------------------------------------------------------------------------------------------------------------------------

void FHEDiskDatabase::genData(uint32_t num_rows, uint32_t num_snp_cols, uint32_t seed)
{
    if (!checkIfDiskDirExists())
    {
        createDiskDir();
    }

    this->num_rows = num_rows;
    this->num_snp_cols = num_snp_cols;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 2);

    this->num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    for (uint32_t i = 0; i < num_snp_cols; i++)
    {
        uint32_t rows_set = 0;

        createColumnDir(i);

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

            saveCtxt(ctxt, i, j);
        }
    }
    this->snp_data_set = true;
    storeDBMetadata();
}

void FHEDiskDatabase::genBinaryPhenoData(uint32_t num_binary_pheno_cols, uint32_t seed)
{
    if (!checkIfDiskDirExists())
    {
        createDiskDir();
    }

    this->num_binary_pheno_cols = num_binary_pheno_cols;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 1);

    for (uint32_t i = 0; i < num_binary_pheno_cols; i++)
    {
        createColumnDirPheno(i, "_binary");

        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = dis(gen);
            }

            helib::Ctxt ctxt = encrypt(ptxt);

            saveCtxtPheno(ctxt, i, j, "_binary");
        }
    }
    this->binary_phenotype_data_set = true;
    storeDBMetadata();
}

void FHEDiskDatabase::genContinuousPhenoData(uint32_t num_continuous_pheno_cols, uint32_t low, uint32_t high, uint32_t seed)
{
    if (!checkIfDiskDirExists())
    {
        createDiskDir();
    }

    this->num_continuous_pheno_cols = num_continuous_pheno_cols;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(low, high);

    for (uint32_t i = 0; i < num_continuous_pheno_cols; i++)
    {
        createColumnDirPheno(i, "_continuous");

        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = dis(gen);
            }

            helib::Ctxt ctxt = encrypt(ptxt);

            saveCtxtPheno(ctxt, i, j, "_continuous");
        }
    }
    this->continuous_phenotype_data_set = true;
    storeDBMetadata();
}

void encryption_thread(FHEDiskDatabase *db, uint32_t start, uint32_t end, vector<vector<int32_t>> &data, uint32_t num_slots, uint32_t num_rows)
{
    for (uint32_t i = start; i < end; i++)
    {
        for (uint32_t j = 0; j < db->num_compressed_rows; j++)
        {
            vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

            uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));

            for (uint32_t k = 0; k < entries_left; k++)
            {
                ptxt[k] = data[i][j * num_slots + k];
            }
            helib::Ctxt ctxt = db->encrypt(ptxt);
            db->saveCtxt(ctxt, i, j);
        }
    }
}

void FHEDiskDatabase::setData(vector<vector<int32_t>> &db)
{
    if (!checkIfDiskDirExists())
    {
        createDiskDir();
    }

    num_snp_cols = db.size();
    if (num_snp_cols == 0)
    {
        throw invalid_argument("ERROR: DB has zero columns! THIS DOES NOT WORK!");
    }

    num_rows = db[0].size();

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    for (uint32_t i = 0; i < num_snp_cols; i++)
    {
        createColumnDir(i);
    }

    vector<thread> threads;

    for (uint32_t i = 0; i < num_threads_encryption; i++)
    {
        uint32_t start = i * num_snp_cols / num_threads_encryption;
        uint32_t end = (i + 1) * num_snp_cols / num_threads_encryption;
        if (i == num_threads_encryption - 1)
        {
            end = num_snp_cols;
        }
        threads.push_back(thread(encryption_thread, this, start, end, ref(db), num_slots, num_rows));
    }

    for (auto &t : threads)
    {
        t.join();
    }

    this->snp_data_set = true;
    storeDBMetadata();
}

void FHEDiskDatabase::multithreadSetData(string vcf_file, uint32_t num_threads)
{
    num_threads_encryption = num_threads;
    FHESIMDDatabase::setData(vcf_file);
}

bool FHEDiskDatabase::checkIfDiskDirExists()
{
    // Check if the disk directory exists
    std::filesystem::path dir(DISK_DIR_FULL);
    return std::filesystem::exists(dir);
}

void FHEDiskDatabase::createDiskDir()
{
    // Create the disk directory
    std::filesystem::path dir(DISK_DIR_FULL);
    std::cout << "Creating disk directory: " << DISK_DIR_FULL << std::endl;
    try
    {
        if (!std::filesystem::exists(dir))
        {
            std::filesystem::create_directories(dir);
        }
        else
        {
            std::cout << "Directory already exists." << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
    }
}

void FHEDiskDatabase::createColumnDir(uint32_t col)
{
    // Create the column directory
    std::filesystem::path dir(this->DISK_DIR_FULL + "/" + std::to_string(col));
    try
    {
        if (!std::filesystem::exists(dir))
        {
            std::filesystem::create_directories(dir);
        }
        else
        {
            std::cout << "Directory already exists. This can cause issues, make sure you delete the database disk before creating a new database" << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
    }
}

void FHEDiskDatabase::createColumnDirPheno(uint32_t col, std::string postfix)
{
    // Create the column directory
    std::filesystem::path dir(this->DISK_DIR_FULL + "/" + std::to_string(col) + postfix);
    try
    {
        if (!std::filesystem::exists(dir))
        {
            std::filesystem::create_directories(dir);
        }
        else
        {
            std::cout << "Directory already exists. This can cause issues, make sure you delete the database disk before creating a new database" << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
    }
}

void FHEDiskDatabase::saveCtxt(helib::Ctxt &ctxt, uint32_t col, uint32_t row)
{
    // Save the ciphertext to disk
    std::string filename = this->DISK_DIR_FULL + "/" + std::to_string(col) + "/" + std::to_string(row) + ".ctxt";
    std::ofstream ofs(filename, std::ios::binary);
    ctxt.writeTo(ofs);
    ofs.close();
}

void FHEDiskDatabase::setGenotype(helib::Ctxt ctxt, uint32_t column, uint32_t compressed_row_index)
{
    saveCtxt(ctxt, column, compressed_row_index);
}

void FHEDiskDatabase::saveCtxtPheno(helib::Ctxt &ctxt, uint32_t col, uint32_t row, string postfix)
{
    // Save the ciphertext to disk
    std::string filename = this->DISK_DIR_FULL + "/" + std::to_string(col) + postfix + "/" + std::to_string(row) + ".ctxt";
    std::ofstream ofs(filename, std::ios::binary);
    ctxt.writeTo(ofs);
    ofs.close();
}

helib::Ctxt FHEDiskDatabase::getGenotype(uint32_t column, uint32_t row) const
{
    if (!snp_data_set)
    {
        std::cout << "SNP data not set." << std::endl;
        exit(1);
    }

    std::string filename = DISK_DIR_FULL + "/" + std::to_string(column) + "/" + std::to_string(row) + ".ctxt";
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs)
    {
        std::cout << "Cannot open file for reading" << std::endl;
        std::cout << "file: " << filename << std::endl;
        exit(1);
    }

    return helib::Ctxt::readFrom(ifs, meta.data->publicKey);
}

helib::Ctxt FHEDiskDatabase::getContinuousPheno(uint32_t column, uint32_t row) const
{
    if (!snp_data_set)
    {
        std::cout << "SNP data not set." << std::endl;
        exit(1);
    }

    std::string filename = DISK_DIR_FULL + "/" + std::to_string(column) + "_continuous/" + std::to_string(row) + ".ctxt";
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs)
    {
        std::cout << "Cannot open file for reading" << std::endl;
        std::cout << "file: " << filename << std::endl;
        exit(1);
    }

    return helib::Ctxt::readFrom(ifs, meta.data->publicKey);
}

helib::Ctxt FHEDiskDatabase::getBinaryPheno(uint32_t column, uint32_t row) const
{
    if (!snp_data_set)
    {
        std::cout << "SNP data not set." << std::endl;
        exit(1);
    }

    std::string filename = DISK_DIR_FULL + "/" + std::to_string(column) + "_binary/" + std::to_string(row) + ".ctxt";
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs)
    {
        std::cout << "Cannot open file for reading" << std::endl;
        std::cout << "file: " << filename << std::endl;
        exit(1);
    }

    return helib::Ctxt::readFrom(ifs, meta.data->publicKey);
}

void FHEDiskDatabase::storeMetadata()
{
    // Save the metadata to disk
    std::string filename = DISK_DIR_FULL + "/" + META_FILEPATH;
    meta.data->saveToFile(filename);
}

void FHEDiskDatabase::storeDBMetadata()
{
    // Save the database metadata to disk
    std::string db_filename = DISK_DIR_FULL + "/" + DB_META_FILE;
    std::ofstream ofs(db_filename);
    ofs << this->num_rows << std::endl;
    ofs << this->num_snp_cols << std::endl;
    ofs << this->num_compressed_rows << std::endl;
    ofs << this->num_continuous_pheno_cols << std::endl;
    ofs << this->num_binary_pheno_cols << std::endl;

    ofs << this->snp_data_set << std::endl;
    ofs << this->continuous_phenotype_data_set << std::endl;
    ofs << this->binary_phenotype_data_set << std::endl;
    ofs.close();
}

void FHEDiskDatabase::loadDBMetadata()
{
    std::string db_filename = DISK_DIR_FULL + "/" + DB_META_FILE;
    std::ifstream ifs(db_filename);
    ifs >> num_rows;
    ifs >> num_snp_cols;
    ifs >> num_compressed_rows;
    ifs >> num_continuous_pheno_cols;
    ifs >> num_binary_pheno_cols;
    ifs >> snp_data_set;
    ifs >> continuous_phenotype_data_set;
    ifs >> binary_phenotype_data_set;
    ifs.close();
}

string FHEDiskDatabase::printDBString(bool with_headers) const
{
    std::stringstream output; // Use a stringstream to build the string

    if (with_headers)
    {
        vector<uint32_t> string_length_count = vector<uint32_t>();

        output << "|";
        for (uint32_t i = 0; i < num_snp_cols; i++)
        {
            output << column_headers[i] << "|";
            string_length_count.push_back(column_headers[i].length());
        }
        output << std::endl;
        output << "--------------";
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {

            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (uint32_t i = 0; i < num_snp_cols; i++)
            {
                temp_storage.push_back(decrypt(getGenotype(i, j)));
            }
            for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
            {

                output << std::endl
                       << "|";

                for (uint32_t i = 0; i < num_snp_cols; i++)
                {
                    if (i > 0)
                    {
                        for (uint32_t space = 0; space < string_length_count[i]; space++)
                        {
                            output << " ";
                        }
                    }

                    output << temp_storage[i][jj];

                    if (i == num_snp_cols - 1)
                    {
                        for (uint32_t space = 0; space < string_length_count[i]; space++)
                        {
                            output << " ";
                        }
                        output << "|";
                    }
                }
            }
        }
    }
    else
    {
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {

            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (uint32_t i = 0; i < num_snp_cols; i++)
            {
                temp_storage.push_back(decrypt(getGenotype(i, j)));
            }
            for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
            {

                output << std::endl
                       << "|";

                for (uint32_t i = 0; i < num_snp_cols; i++)
                {
                    if (i > 0)
                    {
                        output << " ";
                    }

                    output << temp_storage[i][jj];

                    if (i == num_snp_cols - 1)
                    {
                        output << "|";
                    }
                }
            }
        }
    }
    output << std::endl;

    return output.str(); // Return the built string
}

void FHEDiskDatabase::setVCFData(string vcf_file)
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
            if (delimiter_counter < 9)
            {
                continue;
            }
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

void FHEDiskDatabase::setPhenoData(string pheno_file)
{
    vector<string> headers;
    std::ifstream file(pheno_file);
    if (!file.is_open())
    {
        std::cout << "Error opening file: " << pheno_file << std::endl;
        return;
    }
    string header_line;
    getline(file, header_line);
    std::istringstream iss(header_line);
    std::string token;
    while (std::getline(iss, token, ','))
    {
        headers.push_back(token);
    }

    vector<vector<int32_t>> matrix;
    std::string line;
    while (std::getline(file, line))
    {
        std::vector<int32_t> row;
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, ','))
        {
            if (token == "NA")
            {
                row.push_back(-1);
            }
            else
            {
                row.push_back(std::stoi(token));
            }
        }
        matrix.push_back(row);
    }

    // transpose matrix
    vector<vector<int32_t>> transposed_matrix;
    for (uint32_t i = 0; i < matrix[0].size(); i++)
    {
        vector<int32_t> row;
        for (uint32_t j = 0; j < matrix.size(); j++)
        {
            row.push_back(matrix[j][i]);
        }
        transposed_matrix.push_back(row);
    }

    num_binary_pheno_cols = 0;
    num_continuous_pheno_cols = 0;

    for (uint32_t c = 0; c < matrix[0].size(); c++)
    {
        bool is_binary = true;
        for (uint32_t r = 0; r < matrix.size(); r++)
        {
            if (matrix[r][c] != 0 && matrix[r][c] != 1)
            {
                is_binary = false;
                break;
            }
        }
        if (is_binary)
        {
            binary_pheno_headers.push_back(headers[c]);

            createColumnDirPheno(num_binary_pheno_cols, "_binary");


            uint32_t current_row = 0;
            for (uint32_t j = 0; j < num_compressed_rows; j++)
            {
                vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

                uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
                for (uint32_t k = 0; k < entries_left; k++)
                {
                    ptxt[k] = transposed_matrix[c][current_row];
                    current_row += 1;
                }

                helib::Ctxt ctxt = encrypt(ptxt);

                saveCtxtPheno(ctxt, num_binary_pheno_cols, j, "_binary");
            }
            num_binary_pheno_cols += 1;
        }
        else
        {
            continuous_pheno_headers.push_back(headers[c]);

            createColumnDirPheno(num_continuous_pheno_cols, "_continuous");

            uint32_t current_row = 0;
            for (uint32_t j = 0; j < num_compressed_rows; j++)
            {
                vector<unsigned long> ptxt = vector<unsigned long>(num_slots, 0);

                uint32_t entries_left = min(num_slots, num_rows - (j * num_slots));
                for (uint32_t k = 0; k < entries_left; k++)
                {
                    ptxt[k] = transposed_matrix[c][current_row];
                    current_row += 1;
                }

                helib::Ctxt ctxt = encrypt(ptxt);



                saveCtxtPheno(ctxt, num_continuous_pheno_cols, j, "_continuous");
            }
            num_continuous_pheno_cols += 1;
        }
    }
    continuous_phenotype_data_set = true;
    binary_phenotype_data_set = true;
}