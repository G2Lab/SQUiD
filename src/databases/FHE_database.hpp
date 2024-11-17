#pragma once

#include <vector>
#include <string>
#include "FHE_SIMD_database.hpp"
#include <random>
#include <iostream>
#include <thread>
#include <utility>
#include <helib/helib.h>
#include "tools.hpp"
#include "comparator.hpp"
#include "../globals.hpp"

// Database with no packing
class FHEDatabase : public FHESIMDDatabase {
public:
    FHEDatabase(const Params &_params, bool _with_similarity)
        : FHESIMDDatabase(_params, _with_similarity) {
            num_slots = 1;
        }
    ~FHEDatabase() = default;
};