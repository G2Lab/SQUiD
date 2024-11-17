#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

vector<helib::Ctxt> FHESIMDDatabase::PRSQuery(vector<pair<uint32_t, int32_t>> &prs_params) const
{
    vector<helib::Ctxt> scores;

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt sum(meta.data->publicKey);

        for (pair<uint32_t, int32_t> i : prs_params)
        {
            helib::Ctxt temp = getGenotype(i.first, j);
            temp.multByConstant(NTL::ZZX(i.second));

            sum += temp;
        }
        scores.push_back(sum);
    }
    return scores;
}

void process_iteration_prs(FHESIMDDatabase* db,
                           std::vector<helib::Ctxt> &scores,
                           vector<pair<uint32_t, int32_t>> &prs_params,
                           size_t start_idx,
                           size_t end_idx,
                            std::mutex &scores_mutex
                           )
{
    helib::Ctxt score = db->getGenotype(prs_params[start_idx].first, 0);
    score.multByConstant(NTL::ZZX(prs_params[start_idx].second));

    for (size_t i = start_idx + 1; i < end_idx; i++)
    {
        pair<uint32_t, int32_t> param = prs_params[i];
        helib::Ctxt clone = db->getGenotype(param.first, 0);
        clone.multByConstant(NTL::ZZX(param.second));
        score += clone;
    }  

    std::lock_guard<std::mutex> lock(scores_mutex);
    scores.push_back(score);
}

helib::Ctxt FHESIMDDatabase::PRSQueryP(vector<pair<uint32_t, int32_t>> &prs_params, uint32_t num_threads)
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

        if (start_idx >= end_idx)
        {
            continue;
        }

        threads.emplace_back(process_iteration_prs, this,
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

void process_iteration_prs_pp(FHESIMDDatabase* db,
                           std::vector<helib::Ctxt> &scores,
                           vector<pair<uint32_t, int32_t>> &prs_params,
                           size_t start_idx,
                           size_t end_idx,
                           std::mutex &scores_mutex
                           )
{
    for (size_t i = start_idx; i < end_idx; i++)
    {

        helib::Ctxt sum(db->getMeta().data->publicKey);

        for (size_t j = 0; j < prs_params.size(); j++)
        {
            helib::Ctxt temp = db->getGenotype(prs_params[j].first, i);
            temp.multByConstant(NTL::ZZX(prs_params[j].second));
            sum += temp;
        }
        
        {
            std::unique_lock<std::mutex> lock(scores_mutex);
            scores.push_back(sum);
            lock.unlock();
        } // The lock is released here when the block ends
    }  
}

vector<helib::Ctxt> FHESIMDDatabase::PRSQueryPP(vector<pair<uint32_t, int32_t>> &prs_params, uint32_t num_threads)
{
    uint32_t t = num_threads; // You can set this value based on the number of available cores or your requirements

    if (num_compressed_rows < t)
    {
        t = num_compressed_rows;
    }

    size_t chunk_size = num_compressed_rows / t;

    std::vector<std::thread> threads;
    std::vector<helib::Ctxt> scores;

    scores.reserve(num_compressed_rows);

    std::mutex scores_mutex;

    for (size_t i = 0; i < t; i++)
    {
        size_t start_idx = i * chunk_size;
        size_t end_idx = (i == t - 1) ? num_compressed_rows : (i + 1) * chunk_size;

        if (start_idx >= end_idx)
        {
            continue;
        }

        threads.emplace_back(process_iteration_prs_pp, this,
                             std::ref(scores), std::ref(prs_params),
                             start_idx, end_idx, std::ref(scores_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    return scores;
}
