// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef hash_h
#define hash_h

#include <inttypes.h>
#include "MurmurHash3.h"    

typedef uint32_t hash32_t;
typedef uint64_t hash64_t;

union hash_u
{
    hash32_t hash32;
    hash64_t hash64;
};

inline bool hashLessThan(hash_u hash1, hash_u hash2, bool use64)
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


inline hash_u getHash(const char * seq, int length, uint32_t seed, bool use64)
{
    //for ( int i = 0; i < length; i++ ) { cout << *(seq + i); } cout << endl;
    
#ifdef ARCH_32
    char data[use64 ? 8 : 4];
    MurmurHash3_x86_32(seq, length > 16 ? 16 : length, seed, data);
    if ( use64 )
    {
        MurmurHash3_x86_32(seq + 16, length - 16, seed, data + 4);
    }
#else
    char data[16];
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

#endif
