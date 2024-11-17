#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

// Modify Operations
void FHESIMDDatabase::updateOneValue(uint32_t row, uint32_t col, uint32_t value)
{
    helib::Ptxt<helib::BGV> ptxt(meta.data->context);

    uint32_t compressed_row_index = floor(row / num_slots);
    uint32_t row_index = row - (compressed_row_index * num_slots);

    ptxt[row_index] = value;

    helib::Ctxt ctxt(meta.data->publicKey);

    meta.data->publicKey.Encrypt(ctxt, ptxt);

    snp_data[col][compressed_row_index] += ctxt;
}
void FHESIMDDatabase::updateOneRow(uint32_t row, vector<uint32_t> &vals)
{
    for (uint32_t v = 0; v < vals.size(); v++)
    {
        updateOneValue(row, v, vals[v]);
    }
}
void FHESIMDDatabase::insertOneRow(vector<uint32_t> &vals)
{
    uint32_t new_row = num_rows + 1;
    for (uint32_t v = 0; v < vals.size(); v++)
    {
        updateOneValue(new_row, v, vals[v]);
    }
}
void FHESIMDDatabase::deleteRowAddition(uint32_t row)
{
    for (uint32_t v = 0; v < num_snp_cols; v++)
    {
        updateOneValue(row, v, 0);
    }
}

void FHESIMDDatabase::setGenotype(helib::Ctxt ctxt, uint32_t column, uint32_t compressed_row_index)
{
    snp_data[column][compressed_row_index] = ctxt;
}

void FHESIMDDatabase::deleteRowMultiplication(uint32_t row)
{
    helib::Ptxt<helib::BGV> mask(meta.data->context);
    for (uint32_t i = 0; i < num_slots; i++)
    {
        mask[i] = 1;
    }

    uint32_t compressed_row_index = floor(row / num_slots);
    uint32_t row_index = row - (compressed_row_index * num_slots);

    mask[row_index] = 0;

    for (uint32_t c = 0; c < 1; c++)
    {
        helib::Ctxt clone = getGenotype(c,compressed_row_index);
        clone.multByConstant(mask);
        setGenotype(clone, c, compressed_row_index);
    }

    num_deletes += 1;
}
