#pragma once

#include "tools.hpp"

namespace constants
{
    // DEBUG PARAMETERS

    const int DEBUG = 0;
    const int ALPHA = 1; // Maximum number of delete operations before the similarity query fails
    
    // BGV PARAMETERS
    const Params P131(17293, 131, 1, 431);
    const Params BenchParams(65536, 131071, 1, 880);
    const Params BenchParams2(50001, 100003, 1, 880);
}