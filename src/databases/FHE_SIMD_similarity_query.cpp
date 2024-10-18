#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

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

pair<helib::Ctxt, helib::Ctxt> FHESIMDDatabase::similarityQuery(uint32_t target_column, vector<helib::Ctxt> &d, uint32_t threshold) const
{
    // Compute Normalized Score
    if (!with_comparator || !binary_phenotype_data_set)
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
            helib::Ctxt clone = getGenotype(i, j);
            clone -= d[i];
            clone.square();

            temp.push_back(clone);
        }
        normalized_scores.push_back(temp);
    }
    vector<helib::Ctxt> scores = vector<helib::Ctxt>();

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        scores.push_back(addManySafe(normalized_scores[j], meta.data->publicKey));
    }

    if (constants::DEBUG)
    {
        cout << "After scoring:" << endl;
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            print_vector(decrypt(scores[j]));
        }
    }

    helib::Ptxt<helib::BGV> ptxt_threshold(meta.data->context);
    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_threshold[i] = threshold;
    }

    vector<helib::Ctxt> predicate = vector<helib::Ctxt>();
    predicate.reserve(num_compressed_rows);

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt res(meta.data->publicKey);
        comparator->compare(res, scores[j], ptxt_threshold);
        predicate.push_back(std::move(res));
    }

    if (constants::DEBUG)
    {
        cout << "After thresholding:" << endl;
        for (uint32_t j = 0; j < num_compressed_rows; j++)
        {
            print_vector(decrypt(predicate[j]));
        }
    }

    vector<helib::Ctxt> inverse_target_column = vector<helib::Ctxt>();
    inverse_target_column.reserve(num_compressed_rows);

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt inv = getBinaryPheno(target_column, j);
        addOneMod2(inv);
        inverse_target_column.push_back(std::move(inv));
    }

    maskWithNumRows(inverse_target_column);
    maskWithNumRows(predicate);

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        inverse_target_column[j].multiplyBy(predicate[j]);
        predicate[j].multiplyBy(getBinaryPheno(target_column, j));
    }

    helib::Ctxt count_with = addManySafe(predicate, meta.data->publicKey);
    helib::Ctxt count_without = addManySafe(inverse_target_column, meta.data->publicKey);

    count_with = squashCtxtLogTime(count_with);
    count_without = squashCtxtLogTime(count_without);

    return pair(count_with, count_without);
}

void process_iteration_similarity(FHESIMDDatabase *db,
                                  std::vector<helib::Ctxt> &d,
                                  std::vector<helib::Ctxt> &scores,
                                  size_t start_idx,
                                  size_t end_idx,
                                  std::mutex &scores_mutex)
{

    helib::Ctxt score = db->getGenotype(start_idx, 0);
    score -= d[start_idx];
    score.square();
    for (size_t i = start_idx + 1; i < end_idx; i++)
    {
        helib::Ctxt clone = db->getGenotype(i, 0);
        clone -= d[i];
        clone.square();
        score += clone;
    }

    std::lock_guard<std::mutex> lock(scores_mutex);
    scores.push_back(std::move(score));
}

pair<helib::Ctxt, helib::Ctxt> FHESIMDDatabase::similarityQueryP(uint32_t target_column, std::vector<helib::Ctxt> &d, uint32_t threshold, uint32_t num_threads)
{
    if (!with_comparator)
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

        if (start_idx >= end_idx)
        {
            continue;
        }

        threads.emplace_back(process_iteration_similarity, this,
                             std::ref(d), std::ref(scores),
                             start_idx, end_idx, std::ref(scores_mutex));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    helib::Ctxt scores_all = addManySafe(scores, meta.data->publicKey);

    helib::Ptxt<helib::BGV> ptxt_threshold(meta.data->context);
    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_threshold[i] = threshold;
    }
    helib::Ctxt predicate(meta.data->publicKey);
    comparator->compare(predicate, scores_all, ptxt_threshold);

    helib::Ptxt<helib::BGV> mask(meta.data->context);
    for (size_t i = 0; i < num_rows % num_slots; i++)
    {
        mask[i] = 1;
    }
    predicate.multByConstant(mask);

    helib::Ctxt inverse_target = getBinaryPheno(target_column, 0);
    addOneMod2(inverse_target);
    inverse_target.multByConstant(mask);

    inverse_target.multiplyBy(predicate);
    predicate.multiplyBy(getBinaryPheno(target_column, 0));

    helib::Ctxt count_with = predicate;
    helib::Ctxt count_without = inverse_target;

    count_with = squashCtxtLogTime(count_with);
    count_without = squashCtxtLogTime(count_without);

    return pair(count_with, count_without);
}