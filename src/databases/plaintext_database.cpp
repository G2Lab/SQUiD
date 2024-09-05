#include "plaintext_database.hpp"
#include <stdexcept>

using namespace std;

void PlaintextDatabase::genData(uint32_t _num_rows, uint32_t _num_snp_cols, uint32_t seed)
{
    snp_data.resize(_num_snp_cols, vector<int32_t>(_num_rows, 0));
    num_rows = _num_rows;
    num_snp_cols = _num_snp_cols;
    num_compressed_rows = num_rows;

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 2);

    for (size_t i = 0; i < snp_data.size(); ++i)
    {
        for (size_t j = 0; j < snp_data[i].size(); ++j)
        {
            snp_data[i][j] = dis(gen);
        }
    }

    snp_data_set = true;
}

void PlaintextDatabase::genContinuousPhenoData(uint32_t num_continuous_pheno_cols, uint32_t low, uint32_t high, uint32_t seed)
{
    this->num_continuous_pheno_cols = num_continuous_pheno_cols;

    continuous_phenotype_data.resize(num_continuous_pheno_cols, vector<int32_t>(num_rows, 0));

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(low, high);

    for (size_t i = 0; i < continuous_phenotype_data.size(); ++i)
    {
        for (size_t j = 0; j < continuous_phenotype_data[i].size(); ++j)
        {
            continuous_phenotype_data[i][j] = dis(gen);
        }
    }

    continuous_phenotype_data_set = true;
}

void PlaintextDatabase::genBinaryPhenoData(uint32_t num_binary_pheno_cols, uint32_t seed)
{
    this->num_binary_pheno_cols = num_binary_pheno_cols;

    binary_phenotype_data.resize(num_binary_pheno_cols, vector<int32_t>(num_rows, 0));

    random_device rd;
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 1);

    for (size_t i = 0; i < binary_phenotype_data.size(); ++i)
    {
        for (size_t j = 0; j < binary_phenotype_data[i].size(); ++j)
        {
            binary_phenotype_data[i][j] = dis(gen);
        }
    }

    binary_phenotype_data_set = true;
}

void PlaintextDatabase::setData(vector<vector<int32_t>> &db)
{
    snp_data = db;
    num_rows = db[0].size();
    num_snp_cols = db.size();
    num_compressed_rows = num_rows;

    snp_data_set = true;
}

void PlaintextDatabase::setBinaryPhenoData(vector<vector<int32_t>> &db)
{
    binary_phenotype_data = db;
    num_binary_pheno_cols = db.size();

    binary_phenotype_data_set = true;
}

void PlaintextDatabase::setData(string vcf_file)
{
    throw runtime_error("Not implemented");
}

int32_t PlaintextDatabase::countQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const
{
    int32_t count = 0;
    for (size_t i = 0; i < num_rows; i++)
    {
        bool match = true;
        if (conjunctive)
        {
            for (const auto &q : query)
            {
                if (snp_data[q.first][i] != q.second)
                {
                    match = false;
                    break;
                }
            }
        }
        else
        {
            match = false;
            for (const auto &q : query)
            {
                if (snp_data[q.first][i] == q.second)
                {
                    match = true;
                    break;
                }
            }
        }

        if (match)
        {
            count++;
        }
    }
    return count;
}

int32_t PlaintextDatabase::MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const
{
    int32_t count = 0;
    int32_t allele_count = 0;
    for (size_t i = 0; i < num_rows; i++)
    {
        bool match = true;
        if (conjunctive)
        {
            for (const auto &q : query)
            {
                if (snp_data[q.first][i] != q.second)
                {
                    match = false;
                    break;
                }
            }
        }
        else
        {
            match = false;
            for (const auto &q : query)
            {
                if (snp_data[q.first][i] == q.second)
                {
                    match = true;
                    break;
                }
            }
        }

        if (match)
        {
            count += 2;
            allele_count += snp_data[snp][i];
        }
    }
    return count + allele_count * (2 * num_rows); // Hack to return two values in one output
}

vector<int32_t> PlaintextDatabase::PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params) const
{
    vector<int32_t> scores(num_rows, 0);

    for (size_t i = 0; i < num_rows; i++)
    {
        for (const auto &prs : prs_params)
        {
            scores[i] += snp_data[prs.first][i] * prs.second;
        }
    }

    return scores;
}

pair<int32_t, int32_t> PlaintextDatabase::similarityQuery(uint32_t target_column, vector<int32_t> &d, uint32_t threshold) const
{
    int32_t count_with = 0;
    int32_t count_without = 0;

    for (size_t i = 0; i < num_rows; i++)
    {
        uint32_t score = 0;

        for (size_t j = 0; j < d.size(); j++)
        {
            int32_t squared_difference = (snp_data[j][i] - d[j]) * (snp_data[j][i] - d[j]);
            score += squared_difference;
        }

        if (score < threshold)
        {
            if (getBinaryPheno(target_column, i) == 1)
            {
                count_with++;
            }
            else
            {
                count_without++;
            }
        }
    }
    return make_pair(count_with, count_without);
}

int32_t PlaintextDatabase::countingRangeQuery(uint32_t lower, uint32_t upper, uint32_t phenotype_index)
{
    int32_t count = 0;

    for (size_t i = 0; i < num_rows; i++)
    {
        if (getContinuousPheno(phenotype_index, i) >= lower && getContinuousPheno(phenotype_index, i) < upper)
        {
            count++;
        }
    }
    return count;
}

pair<int32_t, int32_t> PlaintextDatabase::MAFRangeQuery(uint32_t snp, uint32_t lower, uint32_t upper, uint32_t phenotype_index)
{
    int32_t count_patients_pass_filter = 0;
    int32_t count_alleles = 0;

    for (size_t i = 0; i < num_rows; i++)
    {
        if (getContinuousPheno(phenotype_index, i) >= lower && getContinuousPheno(phenotype_index, i) < upper)
        {
            count_patients_pass_filter++;
            count_alleles += snp_data[snp][i];
        }
    }
    return make_pair(count_patients_pass_filter, count_alleles);
}

void PlaintextDatabase::printData(bool with_headers) const
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

int32_t PlaintextDatabase::getGenotype(uint32_t i, uint32_t j) const
{
    return snp_data[i][j];
}
int32_t PlaintextDatabase::getContinuousPheno(uint32_t i, uint32_t j) const
{
    return continuous_phenotype_data[i][j];
}
int32_t PlaintextDatabase::getBinaryPheno(uint32_t i, uint32_t j) const
{
    return binary_phenotype_data[i][j];
}

void PlaintextDatabase::printMeta() const
{
    cout << "Number of rows: " << num_rows << endl;
    cout << "Number of snp columns: " << num_snp_cols << endl;
    cout << "Number of binary pheno columns: " << num_binary_pheno_cols << endl;
    cout << "Number of continuous pheno columns: " << num_continuous_pheno_cols << endl;
    cout << "Number of compressed rows: " << num_compressed_rows << endl;
}

void PlaintextDatabase::setGenotype(int32_t data, uint32_t column, uint32_t compressed_row_index)
{
    snp_data[column][compressed_row_index] = data;
}