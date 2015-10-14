// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "hash.h"
#include "MurmurHash3.h"

hash_u getHash(const char * seq, int length)
{
    //for ( int i = 0; i < length; i++ ) { cout << *(seq + i); } cout << endl;
    
    bool use64 = length > 16;
    char data[use64 ? 8 : 4];
    
#ifdef ARCH_32
    MurmurHash3_x86_32(seq, length > 16 ? 16 : length, seed, data);
    if ( use64 )
    {
        MurmurHash3_x86_32(seq + 16, length - 16, seed, data + 4);
    }
#else
    MurmurHash3_x64_128(seq, length, seed, data);
#endif    
    
    hash_u hash;
    
    if ( use64 )
    {
        hash.hash64 = *((hash64_t *)data);
    }
    else
    {
        hash.hash32 = *((hash32_t *)data);
    }
    
    return hash;
}

bool hashLessThan(hash_u hash1, hash_u hash2, bool use64)
{
    if ( use64 )
    {
        return hash1.hash64 < hash2.hash64;
    }
    else
    {
        return hash1.hash32 < hash2.hash32;
    }
}
