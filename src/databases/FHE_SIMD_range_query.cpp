#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

helib::Ctxt FHESIMDDatabase::countingRangeQuery(uint32_t lower, uint32_t upper, uint32_t phenotype_index)
{
    if (!continuous_phenotype_data_set){
        throw invalid_argument("Phenotype data not set");
    }
    if (!with_comparator){
        throw invalid_argument("Count range query not supported without comparator");
    }

    helib::Ptxt<helib::BGV> ptxt_lower(meta.data->context);
    helib::Ptxt<helib::BGV> ptxt_upper(meta.data->context);

    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_lower[i] = lower;
        ptxt_upper[i] = upper;
    }

    vector<helib::Ctxt> predicates = vector<helib::Ctxt>();
    predicates.reserve(num_compressed_rows);

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt lower_predicate(meta.data->publicKey);
        helib::Ctxt upper_predicate(meta.data->publicKey);

        helib::Ctxt phenotype = getContinuousPheno(phenotype_index,j);

        comparator->compare(lower_predicate, phenotype, ptxt_lower);
        comparator->compare(upper_predicate, phenotype, ptxt_upper);

        addOneMod2(lower_predicate);

        upper_predicate *= lower_predicate;
        upper_predicate.cleanUp();

        predicates.push_back(std::move(upper_predicate));
    }

    helib::Ctxt result = addManySafe(predicates, meta.data->publicKey);
    result = squashCtxtLogTime(result);
    return result;
}

pair<helib::Ctxt, helib::Ctxt> FHESIMDDatabase::MAFRangeQuery(uint32_t snp, uint32_t lower, uint32_t upper, uint32_t phenotype_index)
{
    if (!continuous_phenotype_data_set){
        throw invalid_argument("Phenotype data not set");
    }
    if (!with_comparator){
        throw invalid_argument("Count range query not supported without comparator");
    }

    helib::Ptxt<helib::BGV> ptxt_lower(meta.data->context);
    helib::Ptxt<helib::BGV> ptxt_upper(meta.data->context);

    for (uint32_t i = 0; i < num_slots; i++)
    {
        ptxt_lower[i] = lower;
        ptxt_upper[i] = upper;
    }

    vector<helib::Ctxt> predicates = vector<helib::Ctxt>();
    predicates.reserve(num_compressed_rows);

    for (uint32_t j = 0; j < num_compressed_rows; j++)
    {
        helib::Ctxt lower_predicate(meta.data->publicKey);
        helib::Ctxt upper_predicate(meta.data->publicKey);

        helib::Ctxt phenotype = getContinuousPheno(phenotype_index,j);

        comparator->compare(lower_predicate, phenotype, ptxt_lower);
        comparator->compare(upper_predicate, phenotype, ptxt_upper);

        addOneMod2(lower_predicate);

        upper_predicate *= lower_predicate;
        upper_predicate.cleanUp();

        predicates.push_back(upper_predicate);
    }

    helib::Ctxt result = addManySafe(predicates, meta.data->publicKey);

    vector<helib::Ctxt> indv_MAF;
    indv_MAF.reserve(num_compressed_rows);

    for (uint32_t i = 0; i < num_compressed_rows; i++)
    {
        helib::Ctxt clone = getGenotype(snp,i);
        clone *= predicates[i];
        clone.cleanUp();
        
        // Use move semantics to avoid copying the Ctxt object
        indv_MAF.push_back(std::move(clone));
    }

    helib::Ctxt freq = addManySafe(indv_MAF, meta.data->publicKey);
    freq = squashCtxtLogTime(freq);
    result = squashCtxtLogTime(result);

    return pair(freq, result);
}
