#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

void process_iteration_filter(
                              std::vector<helib::Ctxt> &predicates,
                              vector<pair<uint32_t, uint32_t>> &query,
                              FHESIMDDatabase *server_instance,
                              size_t start_idx,
                              size_t end_idx,
                              std::mutex &predicates_mutex);

helib::Ctxt FHESIMDDatabase::MAFQuery(uint32_t snp, bool conjunctive, vector<pair<uint32_t, uint32_t>> &query) const
{
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
    maskWithNumRows(filter_results);

    vector<helib::Ctxt> indv_MAF = vector<helib::Ctxt>();

    for (uint32_t i = 0; i < num_compressed_rows; i++)
    {
        helib::Ctxt clone = getGenotype(snp, i);
        clone *= filter_results[i];
        indv_MAF.push_back(clone);
    }

    helib::Ctxt freq = addManySafe(indv_MAF, meta.data->publicKey);
    helib::Ctxt number_of_patients = addManySafe(filter_results, meta.data->publicKey);

#if COMPRESSED
    freq = squashCtxtWithMask(freq, 0);
    number_of_patients = squashCtxtWithMask(number_of_patients, freq.nAggregates);
#endif

    number_of_patients.multByConstant(NTL::ZZX(2));

    freq += number_of_patients;

    return freq;
}

helib::Ctxt FHESIMDDatabase::MAFQueryP(uint32_t snp, vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads)
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
    helib::Ctxt freq = getGenotype(snp, 0);
    freq *= predicate;
    freq = squashCtxtWithMask(freq, 0);
    helib::Ctxt number_of_patients = squashCtxtWithMask(predicate, freq.nAggregates);
    number_of_patients.multByConstant(NTL::ZZX(2));
    freq += number_of_patients;

    return freq;
}

void process_iteration_MAF_PP(
                              std::vector<helib::Ctxt> &counts,
                              std::vector<helib::Ctxt> &alleles,
                              vector<pair<uint32_t, uint32_t>> &query,
                              uint32_t snp,
                              FHESIMDDatabase *server_instance,
                              size_t start_idx,
                              size_t end_idx,
                              std::mutex &counts_mutex)
{
    helib::Ctxt count(server_instance->getMeta().data->publicKey);
    helib::Ctxt freq(server_instance->getMeta().data->publicKey);
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

        helib::Ctxt target_snp = server_instance->getGenotype(snp, i);

        target_snp *= temp;
        freq += target_snp;
    }

    std::lock_guard<std::mutex> lock(counts_mutex);
    counts.push_back(count);
    alleles.push_back(freq);
}

helib::Ctxt FHESIMDDatabase::MAFQueryPP(uint32_t snp, vector<pair<uint32_t, uint32_t>> &query, uint32_t num_threads)
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
    vector<helib::Ctxt> alleles = vector<helib::Ctxt>();
    counts.reserve(t);
    alleles.reserve(t);

    for (size_t i = 0; i < t; i++)
    {
        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? num_compressed_rows : (i + 1) * chunk_size;

        if (start_idx >= end_idx)
        {
            continue;
        }

        threads.emplace_back(process_iteration_MAF_PP,
                             std::ref(counts), std::ref(alleles), std::ref(query), snp, this,
                             start_idx, end_idx, std::ref(counts_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    helib::Ctxt freq = addManySafe(alleles, meta.data->publicKey);
    helib::Ctxt number_of_patients = addManySafe(counts, meta.data->publicKey); 

    number_of_patients.multByConstant(NTL::ZZX(2));

    freq = squashCtxtWithMask(freq, 0);
    number_of_patients = squashCtxtWithMask(number_of_patients, freq.nAggregates);

    freq += number_of_patients;
    return freq;
}

