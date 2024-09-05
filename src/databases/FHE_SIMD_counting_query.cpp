#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

helib::Ctxt FHESIMDDatabase::countQuery(bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const
{
    if (!snp_data_set)
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
            helib::Ctxt temp = multiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else
    {
        for (uint32_t i = 0; i < num_compressed_rows; i++)
        {
            for (uint32_t j = 0; j < num_columns; j++)
            {
                addOneMod2(cols[i][j]);
            }
        }
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            helib::Ctxt temp = multiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            addOneMod2(filter_results[j]);
        }
    }

    if (constants::DEBUG)
    {
        print_vector(decrypt(filter_results[0]));
    }
    maskWithNumRows(filter_results);
    helib::Ctxt result = addManySafe(filter_results, meta.data->publicKey);

#if COMPRESSED
    result = squashCtxtLogTime(result);
#endif

    return result;
}

void process_iteration_filter(
                              std::vector<helib::Ctxt> &predicates,
                              vector<pair<uint32_t, uint32_t>> &query,
                              FHESIMDDatabase *server_instance,
                              size_t start_idx,
                              size_t end_idx,
                              std::mutex &predicates_mutex)
{
    std::vector<helib::Ctxt> equality_vectors;

    for (size_t i = start_idx; i < end_idx; i++)
    {
        pair<uint32_t, uint32_t> column = query[i];
        equality_vectors.push_back(server_instance->EQTest(column.second, server_instance->getGenotype(column.first, 0)));
    }
    helib::Ctxt predicate = multiplyMany(equality_vectors);

    std::lock_guard<std::mutex> lock(predicates_mutex);
    predicates.push_back(predicate);
}

helib::Ctxt FHESIMDDatabase::countQueryP(vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads)
{
    uint32_t t = num_threads; // You can set this value based on the number of available cores or your requirements

    size_t chunk_size = query.size() / t;

    std::vector<std::thread> threads;

    std::mutex predicates_mutex;

    vector<helib::Ctxt> predicates = vector<helib::Ctxt>();
    predicates.reserve(t);

    for (size_t i = 0; i < t; i++)
    {
        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? query.size() : (i + 1) * chunk_size;

        if (start_idx >= end_idx)
        {
            continue;
        }

        threads.emplace_back(process_iteration_filter,
                             std::ref(predicates), std::ref(query), this,
                             start_idx, end_idx, std::ref(predicates_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }
    
    helib::Ctxt predicate = multiplyMany(predicates);
    return squashCtxtLogTime(predicate);
}

void process_iteration_count_PP(
                              std::vector<helib::Ctxt> &counts,
                              vector<pair<uint32_t, uint32_t>> &query,
                              FHESIMDDatabase *server_instance,
                              size_t start_idx,
                              size_t end_idx,
                              std::mutex &counts_mutex)
{
    helib::Ctxt count(server_instance->getMeta().data->publicKey);
    for (size_t i = start_idx; i < end_idx; i++)
    {
        std::vector<helib::Ctxt> equality_vectors;

        for (size_t j = 0; j < query.size(); j++)
        {
            pair<uint32_t, uint32_t> column = query[j];
            equality_vectors.push_back(server_instance->EQTest(column.second, server_instance->getGenotype(column.first, i)));
        }
        helib::Ctxt temp = multiplyMany(equality_vectors);
        count += temp;
    }

    std::lock_guard<std::mutex> lock(counts_mutex);
    counts.push_back(count);
}

helib::Ctxt FHESIMDDatabase::countQueryPP(vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads)
{
    uint32_t t = num_threads; // You can set this value based on the number of available cores or your requirements

    if (num_compressed_rows < t)
    {
        t = num_compressed_rows;
    }
    
    size_t chunk_size = num_compressed_rows / t;

    std::vector<std::thread> threads;

    std::mutex counts_mutex;

    vector<helib::Ctxt> counts = vector<helib::Ctxt>();
    counts.reserve(t);

    for (size_t i = 0; i < t; i++)
    {
        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? num_compressed_rows : (i + 1) * chunk_size;

        if (start_idx >= end_idx)
        {
            continue;
        }

        threads.emplace_back(process_iteration_count_PP,
                             std::ref(counts), std::ref(query), this,
                             start_idx, end_idx, std::ref(counts_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    helib::Ctxt count = addManySafe(counts, meta.data->publicKey);
    return squashCtxtLogTime(count);
}

