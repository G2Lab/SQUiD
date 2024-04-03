#pragma once

#include "tools.hpp"

namespace constants
{
    // DEBUG PARAMETERS

    const int DEBUG = 0;
    const int ALPHA = 1; // Maximum number of delete operations before the similarity query fails
    
    // BGV PARAMETERS
    const Params N15QP881(32768,65537, 1, 881);
    const Params N14QP438(16384,65537,  1, 438);
    const Params N13QP218(8192, 65537, 1, 218);
    const Params N12QP109(/*m=*/4096, /*p=*/65537, /*r=*/1, /*qbits=*/109);
    const Params P131(17293, 131, 1, 431);
    const Params BenchParams(65536, 131071, 1, 880);
    const Params BenchParams2(50001, 100003, 1, 880);
}
