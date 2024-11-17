#pragma once

#include "./databases/tools.hpp"

namespace constants
{
    // DEBUG PARAMETERS

    const int DEBUG = 0;
    const int ALPHA = 5; // Maximum number of delete operations before the similarity query fails depends on the parameters
    
    // BGV PARAMETERS given as (m, p, r, log q, optional d = 1, optional l = 1)
    // where 
    //m : the cyclotomic order of the ring R;
    //p the plaintext modulus of the BGV scheme;
    //log q : the initial ciphertext modulus of the BGV scheme;
    
    // Insecure but fast parameters used for testing
    const Params Test(2048, 17, 1, 2000);
    
    // Small parameters used for tiny databases
    const Params SmallFastComp(17293, 131, 1, 431);
    
    // Small database which can't run the Similarity or Range queries
    const Params SmallNoSim(8192, 65537, 1, 431);

    // Large parameters used for UKBB scale databases
    const Params Large(99928, 99929, 1, 964);
}