#include "server.hpp"
#include "tools.hpp"

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------

//                                                         HELPER FUNCTIONS

// ------------------------------------------------------------------------------------------------------------------------
helib::Ctxt AddMany(vector<helib::Ctxt> &v);
helib::Ctxt AddManySafe(vector<helib::Ctxt> &v, const helib::PubKey &pk);
helib::Ctxt MultiplyMany(vector<helib::Ctxt> &v);

Server::Server(const Params &_params, bool _with_similarity)
{
    meta(_params);

    num_slots = meta.data->ea.size();
    plaintext_modulus = meta.data->context.getP();

    one_over_two = get_inverse(1, 2, plaintext_modulus);
    neg_three_over_two = get_inverse(-3, 2, plaintext_modulus);
    neg_one_over_two = get_inverse(-1, 2, plaintext_modulus);

    db_set = false;

    with_similarity = _with_similarity;

    // Can take a while to generate polynomials for comparator so we have option to skip it
    if (_with_similarity)
    {
        comparator = unique_ptr<he_cmp::Comparator>(new he_cmp::Comparator(meta.data->context, he_cmp::UNI, 1, 1, meta.data->secretKey, false));
    }
}

void Server::GenData(uint32_t _num_rows, uint32_t _num_cols)
{
    num_rows = _num_rows;
    num_cols = _num_cols;

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    encrypted_db = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_cols; i++)
    {
        vector<helib::Ctxt> cipher_vector = vector<helib::Ctxt>();
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ctxt ctxt = Encrypt(0);

            cipher_vector.push_back(ctxt);
        }
        encrypted_db.push_back(cipher_vector);
    }

    db_set = true;
}

void Server::GenContinuousData(uint32_t _num_rows, uint32_t _low, uint32_t _high)
{
    num_rows = _num_rows;
    num_cols = 1;

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    continuous_db = vector<helib::Ctxt>();
    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ptxt<helib::BGV> ptxt(meta.data->context);

        helib::Ctxt ctxt(meta.data->publicKey);

        for (uint32_t i = 0; i < num_slots; i++)
        {
            ptxt[i] = rand() % (_high - _low + 1) + _low;
        }

        meta.data->publicKey.Encrypt(ctxt, ptxt.getPolyRepr());

        continuous_db.push_back(ctxt);
    }
}

void Server::GenDataDummy(uint32_t _num_rows, uint32_t _num_cols)
{
    num_rows = _num_rows;
    num_cols = _num_cols;

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    encrypted_db = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_cols; i++)
    {
        vector<helib::Ctxt> cipher_vector = vector<helib::Ctxt>();
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ptxt<helib::BGV> ptxt(meta.data->context);

            helib::Ctxt ctxt(meta.data->publicKey);

            ctxt.DummyEncrypt(ptxt.getPolyRepr());

            cipher_vector.push_back(ctxt);
        }
        encrypted_db.push_back(cipher_vector);
    }

    db_set = true;
}

void Server::SetData(vector<vector<uint32_t>> &db)
{
    num_cols = db.size();
    if (num_cols == 0)
    {
        throw invalid_argument("ERROR: DB has zero columns! THIS DOES NOT WORK!");
    }

    num_rows = db[0].size();

    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;

    encrypted_db = vector<vector<helib::Ctxt>>();
    for (uint32_t i = 0; i < num_cols; i++)
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

            helib::Ctxt ctxt = Encrypt(ptxt);

            cipher_vector.push_back(ctxt);
        }
        encrypted_db.push_back(cipher_vector);
    }

    db_set = true;
}

void Server::SetData(string vcf_file)
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

    std::vector<std::vector<uint32_t>> matrix;
    std::string line;
    while (std::getline(file, line))
    {
        // Skip header lines starting with "#"
        if (line[0] == '#')
            continue;

        std::vector<uint32_t> row;

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
            if (token == "./." || token == ".|.")
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
    num_cols = col_counter;

    SetData(matrix);
    file.close();
}

void Server::SetColumnHeaders(vector<string> &headers)
{
    column_headers = vector<string>();
    for (uint32_t i = 0; i < headers.size(); i++)
    {
        column_headers.push_back(headers[i]);
    }
}

// Modify Operations
void Server::UpdateOneValue(uint32_t row, uint32_t col, uint32_t value)
{
    helib::Ptxt<helib::BGV> ptxt(meta.data->context);

    uint32_t compressed_row_index = floor(row / num_slots);
    uint32_t row_index = row - (compressed_row_index * num_slots);

    ptxt[row_index] = value;

    helib::Ctxt ctxt(meta.data->publicKey);

    meta.data->publicKey.Encrypt(ctxt, ptxt);

    encrypted_db[col][compressed_row_index] += ctxt;
}
void Server::UpdateOneRow(uint32_t row, vector<uint32_t> &vals)
{
    for (uint32_t v = 0; v < vals.size(); v++)
    {
        UpdateOneValue(row, v, vals[v]);
    }
}
void Server::InsertOneRow(vector<uint32_t> &vals)
{
    uint32_t new_row = num_rows + 1;
    for (uint32_t v = 0; v < vals.size(); v++)
    {
        UpdateOneValue(new_row, v, vals[v]);
    }
}
void Server::DeleteRowAddition(uint32_t row)
{
    for (uint32_t v = 0; v < num_cols; v++)
    {
        UpdateOneValue(row, v, 0);
    }
}
void Server::DeleteRowMultiplication(uint32_t row)
{
    helib::Ptxt<helib::BGV> mask(meta.data->context);
    for (uint32_t i = 0; i < num_slots; i++)
    {
        mask[i] = 1;
    }

    uint32_t compressed_row_index = floor(row / num_slots);
    uint32_t row_index = row - (compressed_row_index * num_slots);

    mask[row_index] = 0;

    for (uint32_t c = 0; c < num_cols; c++)
    {
        encrypted_db[c][compressed_row_index].multByConstant(mask);
    }

    num_deletes += 1;
}

helib::Ctxt Server::CountQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query)
{
    if (!db_set)
    {
        throw invalid_argument("ERROR: DB needs to be set to run query");
    }

    vector<vector<helib::Ctxt>> cols = filter(query);

    uint32_t num_columns = cols[0].size();

    vector<helib::Ctxt> filter_results;
    if (conjunctive)
    {
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else
    {
        for (uint32_t i = 0; i < num_compressed_rows; i++)
        {
            for (uint32_t j = 0; j < num_columns; j++)
            {
                AddOneMod2(cols[i][j]);
            }
        }
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            AddOneMod2(filter_results[j]);
        }
    }

    if (constants::DEBUG)
    {
        print_vector(Decrypt(filter_results[0]));
    }
    MaskWithNumRows(filter_results);
    helib::Ctxt result = AddManySafe(filter_results, meta.data->publicKey);
    result = SquashCtxtLogTime(result);
    return result;
}

void process_iteration_filter(std::vector<std::vector<helib::Ctxt>> &encrypted_db,
                              std::vector<helib::Ctxt> &predicates,
                              vector<pair<uint32_t, uint32_t>> &query,
                              Server *server_instance,
                              size_t start_idx,
                              size_t end_idx,
                              std::mutex &predicates_mutex)
{
    std::vector<helib::Ctxt> equality_vectors;

    for (size_t i = start_idx; i < end_idx; i++)
    {
        pair<uint32_t, uint32_t> column = query[i];
        equality_vectors.push_back(server_instance->EQTest(column.second, encrypted_db[column.first][0]));
    }
    helib::Ctxt predicate = MultiplyMany(equality_vectors);

    std::lock_guard<std::mutex> lock(predicates_mutex);
    predicates.push_back(predicate);
}

helib::Ctxt Server::CountQueryP(vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads)
{
    uint32_t t = num_threads; // You can set this value based on the number of available cores or your requirements

    size_t chunk_size = query.size() / t;

    std::vector<std::thread> threads;

    std::mutex predicates_mutex;

    vector<helib::Ctxt> predicates = vector<helib::Ctxt>();
    predicates.reserve(query.size());

    for (size_t i = 0; i < t; i++)
    {
        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? query.size() : (i + 1) * chunk_size;

        threads.emplace_back(process_iteration_filter, std::ref(encrypted_db),
                             std::ref(predicates), std::ref(query), this,
                             start_idx, end_idx, std::ref(predicates_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }
    helib::Ctxt predicate = MultiplyMany(predicates);
    return SquashCtxtLogTime(predicate);
}

helib::Ctxt Server::MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query)
{
    vector<vector<helib::Ctxt>> cols = filter(query);
    uint32_t num_columns = cols[0].size();

    vector<helib::Ctxt> filter_results;
    if (conjunctive)
    {
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else
    {
        for (uint32_t i = 0; i < num_compressed_rows; i++)
        {
            for (uint32_t j = 0; j < num_columns; j++)
            {
                AddOneMod2(cols[i][j]);
            }
        }
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            AddOneMod2(filter_results[j]);
        }
    }
    MaskWithNumRows(filter_results);

    vector<helib::Ctxt> indv_MAF = vector<helib::Ctxt>();

    for (uint32_t i = 0; i < num_compressed_rows; i++)
    {
        helib::Ctxt clone = encrypted_db[snp][i];
        clone *= filter_results[i];
        indv_MAF.push_back(clone);
    }

    helib::Ctxt freq = AddManySafe(indv_MAF, meta.data->publicKey);
    helib::Ctxt number_of_patients = AddManySafe(filter_results, meta.data->publicKey);

    freq = SquashCtxtWithMask(freq, 0);
    number_of_patients = SquashCtxtWithMask(number_of_patients, 1);

    number_of_patients.multByConstant(NTL::ZZX(2));

    freq += number_of_patients;

    return freq;
}

helib::Ctxt Server::MAFQueryP(uint32_t snp, vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads)
{
    uint32_t t = num_threads; // You can set this value based on the number of available cores or your requirements

    size_t chunk_size = query.size() / t;

    std::vector<std::thread> threads;

    std::mutex predicates_mutex;

    vector<helib::Ctxt> predicates = vector<helib::Ctxt>();
    predicates.reserve(query.size());

    for (size_t i = 0; i < t; i++)
    {
        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? query.size() : (i + 1) * chunk_size;

        threads.emplace_back(process_iteration_filter, std::ref(encrypted_db),
                             std::ref(predicates), std::ref(query), this,
                             start_idx, end_idx, std::ref(predicates_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }
    helib::Ctxt predicate = MultiplyMany(predicates);

    helib::Ctxt freq = encrypted_db[snp][0];
    freq *= predicate;

    freq = SquashCtxtWithMask(freq, 0);
    helib::Ctxt number_of_patients = SquashCtxtWithMask(predicate, 1);

    number_of_patients.multByConstant(NTL::ZZX(2));


    freq += number_of_patients;
    return freq;
}

vector<helib::Ctxt> Server::PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params)
{
    vector<helib::Ctxt> scores;

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt sum(meta.data->publicKey);

        for (pair<uint32_t, int32_t> i : prs_params)
        {
            helib::Ctxt temp = encrypted_db[i.first][j];
            temp.multByConstant(NTL::ZZX(i.second));

            sum += temp;
        }
        scores.push_back(sum);
    }
    return scores;
}

void process_iteration_prs(std::vector<std::vector<helib::Ctxt>> &encrypted_db,
                           std::vector<helib::Ctxt> &scores,
                           vector<pair<uint32_t, int32_t>> &prs_params,
                           size_t start_idx,
                           size_t end_idx,
                            std::mutex &scores_mutex
                           )
{
    helib::Ctxt score = encrypted_db[prs_params[start_idx].first][0];
    score.multByConstant(NTL::ZZX(prs_params[start_idx].second));

    for (size_t i = start_idx + 1; i < end_idx; i++)
    {
        pair<uint32_t, int32_t> param = prs_params[i];
        helib::Ctxt clone = encrypted_db[param.first][0];
        clone.multByConstant(NTL::ZZX(param.second));
        score += clone;
    }  

    std::lock_guard<std::mutex> lock(scores_mutex);
    scores.push_back(score);
}

helib::Ctxt Server::PRSQueryP(vector<pair<uint32_t, int32_t>> &prs_params, uint32_t num_threads)
{
    uint32_t t = num_threads; // You can set this value based on the number of available cores or your requirements

    size_t chunk_size = prs_params.size() / t;

    std::vector<std::thread> threads;

    vector<helib::Ctxt> scores = vector<helib::Ctxt>();
    scores.reserve(t);

    std::mutex scores_mutex;


    for (size_t i = 0; i < t; i++)
    {
        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? prs_params.size() : (i + 1) * chunk_size;

        threads.emplace_back(process_iteration_prs, std::ref(encrypted_db),
                             std::ref(scores), std::ref(prs_params),
                             start_idx, end_idx, std::ref(scores_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    helib::Ctxt scores_all = scores[0];
    for (size_t i = 1; i < t; i++)
    {
        scores_all += scores[i];
    }
    return scores_all;
}

pair<helib::Ctxt, helib::Ctxt> Server::SimilarityQuery(uint32_t target_column, vector<helib::Ctxt> &d, uint32_t threshold)
{
    // Compute Normalized Score
    if (!with_similarity)
    {
        std::cout << "Server not setup to run similarity queries" << std::endl;
        throw "Invalid setup";
    }

    if (num_deletes > constants::ALPHA)
    {
        std::cout << "Too many deletes have been performed. Cannot run similarity query" << std::endl;
        std::cout << "The data owner needs to refresh the ciphertexts" << std::endl;
        throw "Too many deletes";
    }

    vector<vector<helib::Ctxt>> normalized_scores = vector<vector<helib::Ctxt>>();

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        vector<helib::Ctxt> temp = vector<helib::Ctxt>();
        for (size_t i = 0; i < d.size(); i++)
        {
            helib::Ctxt clone = encrypted_db[i][j];
            clone -= d[i];
            clone.square();
            
            clone = clone.cleanUp();

            temp.push_back(clone);
        }
        normalized_scores.push_back(temp);
    }
    vector<helib::Ctxt> scores = vector<helib::Ctxt>();

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        scores.push_back(AddManySafe(normalized_scores[j], meta.data->publicKey));
    }
    if (constants::DEBUG)
    {
        cout << "After scoring:" << endl;
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            print_vector(Decrypt(scores[j]));
        }
    }

    helib::Ptxt<helib::BGV> ptxt_threshold(meta.data->context);
    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_threshold[i] = threshold;
    }

    vector<helib::Ctxt> predicate = vector<helib::Ctxt>();
    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt res(meta.data->publicKey);
        comparator->compare(res, scores[j], ptxt_threshold);
        predicate.push_back(res);
    }

    if (constants::DEBUG)
    {
        cout << "After thresholding:" << endl;
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            print_vector(Decrypt(predicate[j]));
        }
    }

    vector<helib::Ctxt> inverse_target_column = vector<helib::Ctxt>();
    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt inv = encrypted_db[target_column][j];
        AddOneMod2(inv);
        inverse_target_column.push_back(inv);
    }

    MaskWithNumRows(inverse_target_column);
    MaskWithNumRows(predicate);

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        inverse_target_column[j].multiplyBy(predicate[j]);
        inverse_target_column[j] = inverse_target_column[j].cleanUp();

        predicate[j].multiplyBy(encrypted_db[target_column][j]);
        predicate[j] = predicate[j].cleanUp();

    }

    helib::Ctxt count_with = AddManySafe(predicate, meta.data->publicKey);
    helib::Ctxt count_without = AddManySafe(inverse_target_column, meta.data->publicKey);

    count_with = SquashCtxtLogTime(count_with);
    count_without = SquashCtxtLogTime(count_without);

    return pair(count_with, count_without);
}

void process_iteration_similarity(std::vector<std::vector<helib::Ctxt>> &encrypted_db,
                                  std::vector<helib::Ctxt> &d,
                                  std::vector<helib::Ctxt> &scores,
                                  size_t start_idx,
                                  size_t end_idx,
                                  std::mutex &scores_mutex)
{

    helib::Ctxt score = encrypted_db[start_idx][0];
    score -= d[start_idx];
    score.square();
    score.cleanUp();
    for (size_t i = start_idx; i < end_idx; i++)
    {
        helib::Ctxt clone = encrypted_db[i][0];
        clone -= d[i];
        clone.square();
        clone.cleanUp();
        score += clone;
    }

    std::lock_guard<std::mutex> lock(scores_mutex);
    scores.push_back(score);
}

pair<helib::Ctxt, helib::Ctxt> Server::SimilarityQueryP(uint32_t target_column, std::vector<helib::Ctxt> &d, uint32_t threshold, uint32_t num_threads)
{
    if (!with_similarity)
    {
        std::cout << "Server not setup to run similarity queries" << std::endl;
        throw "Invalid setup";
    }

    uint32_t t = num_threads;
    uint32_t num_snps = d.size();

    size_t chunk_size = num_snps / t;

    vector<helib::Ctxt> scores = vector<helib::Ctxt>();
    scores.reserve(t);

    std::mutex scores_mutex;

    std::vector<std::thread> threads;
    for (size_t i = 0; i < t; i++)
    {

        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? num_snps : (i + 1) * chunk_size;

        threads.emplace_back(process_iteration_similarity, std::ref(encrypted_db),
                             std::ref(d), std::ref(scores),
                             start_idx, end_idx, std::ref(scores_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    helib::Ctxt scores_all = AddManySafe(scores, meta.data->publicKey);

    helib::Ptxt<helib::BGV> ptxt_threshold(meta.data->context);
    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_threshold[i] = threshold;
    }

    helib::Ctxt predicate(meta.data->publicKey);
        
    comparator->compare(predicate, scores_all, ptxt_threshold);

    helib::Ctxt inverse_predicate = predicate;
    AddOneMod2(inverse_predicate);

    predicate *= encrypted_db[target_column][0];
    inverse_predicate *= encrypted_db[target_column][0];

    helib::Ctxt count_with = SquashCtxtLogTime(predicate);
    helib::Ctxt count_without = SquashCtxtLogTime(inverse_predicate);

    return pair(count_with, count_without);
}

helib::Ctxt Server::CountingRangeQuery(uint32_t  lower, uint32_t  upper)
{
    helib::Ptxt<helib::BGV> ptxt_lower(meta.data->context);
    helib::Ptxt<helib::BGV> ptxt_upper(meta.data->context);

    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_lower[i] = lower;
        ptxt_upper[i] = upper;
    }

    vector<helib::Ctxt> predicates = vector<helib::Ctxt>();

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt lower_predicate(meta.data->publicKey);
        helib::Ctxt upper_predicate(meta.data->publicKey);

        comparator->compare(lower_predicate, continuous_db[j], ptxt_lower);
        comparator->compare(upper_predicate, continuous_db[j], ptxt_upper);

        AddOneMod2(lower_predicate);

        upper_predicate *= lower_predicate;
        upper_predicate.cleanUp();

        predicates.push_back(upper_predicate);
    }

    helib::Ctxt result = AddManySafe(predicates, meta.data->publicKey);
    result = SquashCtxtLogTime(result);
    return result;
}

pair<helib::Ctxt, helib::Ctxt> Server::MAFRangeQuery(uint32_t  snp, uint32_t  lower, uint32_t  upper)
{
    helib::Ptxt<helib::BGV> ptxt_lower(meta.data->context);
    helib::Ptxt<helib::BGV> ptxt_upper(meta.data->context);

    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_lower[i] = lower;
        ptxt_upper[i] = upper;
    }

    vector<helib::Ctxt> predicates = vector<helib::Ctxt>();

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt lower_predicate(meta.data->publicKey);
        helib::Ctxt upper_predicate(meta.data->publicKey);

        comparator->compare(lower_predicate, continuous_db[j], ptxt_lower);
        comparator->compare(upper_predicate, continuous_db[j], ptxt_upper);

        AddOneMod2(lower_predicate);

        upper_predicate *= lower_predicate;
        upper_predicate.cleanUp();

        predicates.push_back(upper_predicate);
    }

    helib::Ctxt result = AddManySafe(predicates, meta.data->publicKey);
    result = SquashCtxtLogTime(result);

    vector<helib::Ctxt> indv_MAF = vector<helib::Ctxt>();

    for (uint32_t i = 0; i < num_compressed_rows; i++)
    {
        helib::Ctxt clone = encrypted_db[snp][i];
        clone *= predicates[i];
        clone.cleanUp();
        indv_MAF.push_back(clone);
    }

    helib::Ctxt freq = AddManySafe(indv_MAF, meta.data->publicKey);
    freq = SquashCtxtLogTime(freq);
    return pair(freq, result);
}


void Server::AddOneMod2(helib::Ctxt &a)
{
    //   0 -> 1
    //   1 -> 0
    // f(x)-> -x+1

    a.negate();
    a.addConstant(NTL::ZZX(1));
}

helib::Ctxt MultiplyMany(vector<helib::Ctxt> &v)
{
    uint32_t num_entries = v.size();
    uint32_t depth = ceil(log2(num_entries));

    for (uint32_t d = 0; d < depth; d++)
    {
        uint32_t jump_factor = pow(2, d);
        uint32_t skip_factor = 2 * jump_factor;

        for (uint32_t i = 0; i < num_entries; i += skip_factor)
        {
            v[i].multiplyBy(v[i + jump_factor]);
        }
    }
    return v[0];
}

helib::Ctxt AddMany(vector<helib::Ctxt> &v)
{
    uint32_t num_entries = v.size();
    uint32_t depth = ceil(log2(num_entries));

    for (uint32_t d = 0; d < depth; d++)
    {
        uint32_t jump_factor = pow(2, d);
        uint32_t skip_factor = 2 * jump_factor;

        for (uint32_t i = 0; i < num_entries; i += skip_factor)
        {
            v[i] += v[i + jump_factor];
        }
    }
    return v[0];
}

helib::Ctxt AddManySafe(vector<helib::Ctxt> &v, const helib::PubKey &pk)
{
    helib::Ctxt result(pk);

    for (uint32_t i = 0; i < v.size(); i++)
    {
        result += v[i];
    }
    return result;
}

helib::Ctxt Server::SquashCtxt(helib::Ctxt &ciphertext, uint32_t num_data_elements)
{
    const helib::EncryptedArray &ea = meta.data->context.getEA();

    helib::Ctxt result = ciphertext;

    for (uint32_t i = 1; i < num_data_elements; i++)
    {
        ea.rotate(ciphertext, -(1));
        result += ciphertext;
    }
    return result;
}

helib::Ctxt Server::SquashCtxtLogTimePower2(helib::Ctxt &ciphertext)
{
    const helib::EncryptedArray &ea = meta.data->context.getEA();

    uint32_t depth = floor(log2(num_slots));

    for (int d = depth - 1; d >= 0; d--)
    {
        int32_t shift = 1 << d;
        helib::Ctxt clone = ciphertext;
        ea.rotate(clone, (-shift));
        ciphertext += clone;
        cout << "Step " << d << " capacity:" << ciphertext.capacity() << endl;
    }

    return ciphertext;
}

helib::Ctxt Server::SquashCtxtLogTime(helib::Ctxt &ciphertext)
{
    const helib::EncryptedArray &ea = meta.data->context.getEA();

    uint32_t depth = floor(log2(num_slots));

    helib::Ptxt<helib::BGV> mask(meta.data->context);
    helib::Ptxt<helib::BGV> inverse_mask(meta.data->context);

    uint32_t largest_power_of_two_less_than_or_equal_two_slotsize = 1 << depth;
    for (uint32_t i = 0; i < num_slots; i++)
    {
        if (i < largest_power_of_two_less_than_or_equal_two_slotsize)
        {
            mask[i] = 1;
            inverse_mask[i] = 0;
        }
        else
        {
            mask[i] = 0;
            inverse_mask[i] = 1;
        }
    }
    helib::Ctxt far_end = ciphertext;
    far_end.multByConstant(inverse_mask);

    ea.rotate(far_end, -largest_power_of_two_less_than_or_equal_two_slotsize);

    ciphertext.multByConstant(mask);
    ciphertext += far_end;

    for (int d = depth - 1; d >= 0; d--)
    {
        int32_t shift = 1 << d;
        helib::Ctxt clone = ciphertext;
        ea.rotate(clone, (-shift));
        ciphertext += clone;
    }

    return ciphertext;
}

helib::Ctxt Server::SquashCtxtWithMask(helib::Ctxt &ciphertext, uint32_t index)
{
    ciphertext = SquashCtxtLogTime(ciphertext);

    const helib::EncryptedArray &ea = meta.data->context.getEA();
    if (index != 0)
    {
        ea.rotate(ciphertext, index);
    }
    helib::Ptxt<helib::BGV> mask(meta.data->context);
    mask[index] = 1;
    ciphertext.multByConstant(mask);

    return ciphertext;
}

void Server::CtxtExpand(helib::Ctxt &ciphertext)
{
    const helib::EncryptedArray &ea = meta.data->context.getEA();

    uint32_t depth = floor(log2(num_slots));

    helib::Ptxt<helib::BGV> mask(meta.data->context);

    uint32_t largest_power_of_two_less_than_or_equal_two_slotsize = 1 << depth;
    for (uint32_t i = 0; i < num_slots - largest_power_of_two_less_than_or_equal_two_slotsize; i++)
    {
        mask[i] = 1;
    }
    for (uint32_t d = 0; d < depth; d++)
    {
        uint32_t shift = 1 << d;
        helib::Ctxt clone = ciphertext;
        ea.rotate(clone, (shift));
        ciphertext += clone;
    }

    helib::Ctxt clone = ciphertext;
    clone.multByConstant(mask);
    ea.rotate(clone, largest_power_of_two_less_than_or_equal_two_slotsize);
    ciphertext += clone;
}

void Server::MaskWithNumRows(vector<helib::Ctxt> &ciphertexts)
{
    helib::Ptxt<helib::BGV> mask(meta.data->context);
    for (size_t i = 0; i < num_rows % num_slots; i++)
    {
        mask[i] = 1;
    }
    ciphertexts.back().multByConstant(mask);
}

helib::Ctxt Server::EQTest(unsigned long a, helib::Ctxt &b)
{
    helib::Ctxt clone = b;
    helib::Ctxt result = b;

    switch (a)
    {
    case 0:
    {
        // f(x) = x^2 / 2 - 3/2 x + 1
        //  0 -> 1
        //  1 -> 0
        //  2 -> 0

        result.square();

        result.multByConstant(NTL::ZZX(one_over_two));
        clone.multByConstant(NTL::ZZX(neg_three_over_two));

        result += clone;

        result.addConstant(NTL::ZZX(1));

        return result;
    }
    case 1:
    {
        // f(x) = -x^2 + 2x
        //  0 -> 0
        //  1 -> 1
        //  2 -> 0
        clone.square();

        result.multByConstant(NTL::ZZX(2));

        result -= clone;
        return result;
    }
    case 2:
    {
        // f(x) = x^2 / 2 - x / 2
        //  0 -> 0
        //  1 -> 0
        //  2 -> 1

        result.square();

        result.multByConstant(NTL::ZZX(one_over_two));
        clone.multByConstant(NTL::ZZX(neg_one_over_two));

        result += clone;
        return result;
    }
    default:
        cout << "Can't use a value of a other than 0, 1, or 2" << endl;
        throw invalid_argument("ERROR: invalid value for EQTest");
    }
}

vector<vector<helib::Ctxt>> Server::filter(vector<pair<uint32_t, uint32_t>> &query)
{
    vector<vector<helib::Ctxt>> feature_cols;

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        vector<helib::Ctxt> indv_vector;
        for (pair<uint32_t, uint32_t> i : query)
        {
            indv_vector.push_back(EQTest(i.second, encrypted_db[i.first][j]));

            if (constants::DEBUG == 3)
            {
                cout << "checking equality to " << i.second << endl;
                cout << "original:";
                print_vector(Decrypt(encrypted_db[i.first][j]));
                cout << "result  :";
                print_vector(Decrypt(EQTest(i.second, encrypted_db[i.first][j])));
            }
        }
        feature_cols.push_back(indv_vector);
    }
    return feature_cols;
}

vector<long> Server::Decrypt(helib::Ctxt ctxt)
{
    helib::Ptxt<helib::BGV> new_plaintext_result(meta.data->context);
    meta.data->secretKey.Decrypt(new_plaintext_result, ctxt);

    vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();

    vector<long> result = vector<long>(num_slots);

    for (uint32_t i = 0; i < num_slots; i++)
    {
        result[i] = (long)poly_mod_result[i];
    }

    return result;
}

helib::Ptxt<helib::BGV> Server::DecryptPlaintext(helib::Ctxt ctxt)
{

    if (constants::DEBUG && ctxt.capacity() < 2)
    {
        cout << "NOISE BOUNDS EXCEEDED!!!" << endl;
    }

    helib::Ptxt<helib::BGV> new_plaintext_result(meta.data->context);
    meta.data->secretKey.Decrypt(new_plaintext_result, ctxt);

    return new_plaintext_result;
}

helib::Ctxt Server::Encrypt(unsigned long a)
{
    vector<unsigned long> a_vec = vector<unsigned long>();
    for (size_t i = 0; i < num_slots; i++)
    {
        a_vec.push_back(a);
    }

    return Encrypt(a_vec);
}

helib::Ctxt Server::Encrypt(vector<unsigned long> a)
{
    if (a.size() > num_slots)
    {
        throw invalid_argument("Trying to encrypt vector with too many elements");
    }
    helib::Ptxt<helib::BGV> ptxt(meta.data->context);

    for (size_t i = 0; i < a.size(); ++i)
    {
        ptxt[i] = a[i];
    }

    helib::Ctxt ctxt(meta.data->publicKey);

    meta.data->publicKey.Encrypt(ctxt, ptxt);

    return ctxt;
}

helib::Ctxt Server::EncryptSK(unsigned long a)
{
    vector<unsigned long> a_vec = vector<unsigned long>();
    for (size_t i = 0; i < num_slots; i++)
    {
        a_vec.push_back(a);
    }

    return EncryptSK(a_vec);
}

helib::Ctxt Server::EncryptSK(vector<unsigned long> a)
{
    if (a.size() > num_slots)
    {
        throw invalid_argument("Trying to encrypt vector with too many elements");
    }
    helib::Ptxt<helib::BGV> ptxt(meta.data->context);

    for (size_t i = 0; i < a.size(); ++i)
    {
        ptxt[i] = a[i];
    }

    helib::Ctxt ctxt(meta.data->secretKey);

    EncodedPtxt eptxt;
    ptxt.encode(eptxt);

    meta.data->secretKey.Encrypt(ctxt, eptxt);

    return ctxt;
}

helib::Ctxt Server::GetAnyElement()
{
    return encrypted_db[0][0];
}

void Server::PrintContext()
{
    meta.data->context.printout();
    cout << endl;
    cout << "Security: " << meta.data->context.securityLevel() << endl;
    cout << "Num slots: " << num_slots << endl;
    cout << "Num rows:" << num_rows << endl;
    cout << "Num compressed rows: " << num_compressed_rows << endl;
    cout << "Num cols: " << num_cols << endl;
}

void Server::PrintEncryptedDB(bool with_headers)
{
    if (with_headers)
    {
        vector<uint32_t> string_length_count = vector<uint32_t>();

        cout << "|";
        for (uint32_t i = 0; i < num_cols; i++)
        {
            cout << column_headers[i] << "|";
            string_length_count.push_back(column_headers[i].length());
        }
        cout << endl;
        cout << "--------------";
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {

            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (uint32_t i = 0; i < num_cols; i++)
            {
                temp_storage.push_back(Decrypt(encrypted_db[i][j]));
            }
            for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
            {

                cout << endl
                     << "|";

                for (uint32_t i = 0; i < num_cols; i++)
                {
                    if (i > 0)
                    {
                        for (uint32_t space = 0; space < string_length_count[i]; space++)
                        {
                            cout << " ";
                        }
                    }

                    cout << temp_storage[i][jj];

                    if (i == num_cols - 1)
                    {
                        for (uint32_t space = 0; space < string_length_count[i]; space++)
                        {
                            cout << " ";
                        }
                        cout << "|";
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
            for (uint32_t i = 0; i < num_cols; i++)
            {
                temp_storage.push_back(Decrypt(encrypted_db[i][j]));
            }
            for (uint32_t jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++)
            {

                cout << endl
                     << "|";

                for (uint32_t i = 0; i < num_cols; i++)
                {
                    if (i > 0)
                    {
                        cout << " ";
                    }

                    cout << temp_storage[i][jj];

                    if (i == num_cols - 1)
                    {
                        cout << "|";
                    }
                }
            }
        }
    }
    cout << endl;
}

uint32_t Server::GetSlotSize()
{
    return num_slots;
}

uint32_t Server::GetCompressedRows()
{
    return num_compressed_rows;
}

uint32_t Server::GetCols()
{
    return num_cols;
}

vector<string> Server::GetHeaders(){
    return column_headers;
}


// IMPORTED FROM HELIB SOURCE CODE
inline long estimateCtxtSize(const helib::Context &context, long offset)
{
    // Return in bytes.

    // We assume that the size of each element in the DCRT is BINIO_64BIT

    // sizeof(BINIO_EYE_CTXT_BEGIN) = 4;
    // BINIO_32BIT = 4
    // sizeof(long) = BINIO_64BIT = 8
    // xdouble = s * sizeof(long) = 2 * BINIO_64BIT = 16

    // We assume that primeSet after encryption is context.ctxtPrimes
    // We assume we have exactly 2 parts after encryption
    // We assume that the DCRT prime set is the same as the ctxt one

    long size = 0;

    // Header metadata
    size += 24;

    // Begin eye-catcher
    size += 4;

    // Begin Ctxt metadata
    // 64 = header_size = ptxtSpace (long) + intFactor (long) + ptxtMag (xdouble)
    //                    + ratFactor (xdouble) + noiseBound (xdouble)
    size += 64;

    // primeSet.write(str);
    // size of set (long) + each prime (long)
    size += 8 + context.getCtxtPrimes().card() * 8;

    // Begin Ctxt content size
    // write_raw_vector(str, parts);
    // Size of the parts vector (long)
    size += 8;

    long part_size = 0;
    // Begin CtxtPart size

    // skHandle.write(str);
    // powerOfS (long) + powerOfX (long) + secretKeyID (long)
    part_size += 24;

    // Begin DCRT size computation

    // this->DoubleCRT::write(str);
    // map.getIndexSet().write(str);
    // size of set (long) + each prime (long)
    part_size += 8 + context.getCtxtPrimes().card() * 8;

    // DCRT data write as write_ntl_vec_long(str, map[i]);
    // For each prime in the ctxt modulus chain
    //    size of DCRT column (long) + size of each element (long) +
    //    size of all the slots (column in DCRT) (PhiM long elements)
    long dcrt_size = (8 + 8 * context.getPhiM()) * context.getCtxtPrimes().card();

    part_size += dcrt_size;

    // End DCRT size
    // End CtxtPart size

    size += 2 * part_size; // 2 * because we assumed 2 parts
    // End Ctxt content size

    // End eye-catcher
    size += 4;

    return size + offset;
}

uint32_t Server::StorageOfOneElement()
{
    if (!db_set)
    {
        throw invalid_argument("ERROR: DB needs to be set to get storage cost");
    }

    return estimateCtxtSize(meta.data->context, 0);
}

template <typename T, typename Allocator>
void print_vector(const vector<T, Allocator> &vect, int num_entries)
{
    cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++)
    {
        cout << ", " << vect[i];
    }
    cout << endl;
}
