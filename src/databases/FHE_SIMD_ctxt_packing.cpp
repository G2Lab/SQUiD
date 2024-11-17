#include "FHE_SIMD_database.hpp"
#include "tools.hpp"

using namespace std;

helib::Ctxt FHESIMDDatabase::squashCtxt(helib::Ctxt &ciphertext) const
{
    const helib::EncryptedArray &ea = meta.data->context.getEA();
    helib::Ctxt sum_ciphertext = ciphertext;

    long rotation = 1;
    while (rotation < num_slots) {
        helib::Ctxt rotated_ciphertext = ciphertext;
        ea.rotate(rotated_ciphertext, rotation);
        sum_ciphertext += rotated_ciphertext;
        rotation *= 2;
    }

    return sum_ciphertext;
}

helib::Ctxt FHESIMDDatabase::squashCtxtLogTimePower2(helib::Ctxt &ciphertext) const
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

helib::Ctxt FHESIMDDatabase::squashCtxtLogTime(helib::Ctxt &ciphertext) const
{
    const helib::EncryptedArray &ea = meta.data->context.getEA();

    uint32_t depth = floor(log2(num_slots));

    helib::Ptxt<helib::BGV> mask(meta.data->context);
    helib::Ptxt<helib::BGV> inverse_mask(meta.data->context);

    int32_t chunkingFactor = getChunkingFactor();

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

    uint32_t level = 0;

    for (int d = depth - 1; d >= 0; d--)
    {
        int32_t shift = 1 << d;
        helib::Ctxt clone = ciphertext;
        ea.rotate(clone, (-shift));
        ciphertext += clone;

        if (level == chunkingFactor || ciphertext.bitCapacity() < EARLY_TERM_NOISE_THRES)
        {
            ciphertext.nAggregates = shift;
            break;
        }

        level += 1;
    }

    return ciphertext;
}

helib::Ctxt FHESIMDDatabase::squashCtxtWithMask(helib::Ctxt &ciphertext, uint32_t index) const
{
    ciphertext = squashCtxtLogTime(ciphertext);

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

void FHESIMDDatabase::ctxtExpand(helib::Ctxt &ciphertext) const
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

void FHESIMDDatabase::maskWithNumRows(vector<helib::Ctxt> &ciphertexts) const
{
    helib::Ptxt<helib::BGV> mask(meta.data->context);
    for (size_t i = 0; i < num_rows % num_slots; i++)
    {
        mask[i] = 1;
    }
    if (num_rows % num_slots != 0)
    {
        ciphertexts.back().multByConstant(mask);
    }
}