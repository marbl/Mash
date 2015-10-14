// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "HashSet.h"

void HashSet::erase(hash_u hash)
{
    if ( use64 )
    {
        hashes64.erase(hash.hash64);
    }
    else
    {
        hashes32.erase(hash.hash32);
    }
}

void HashSet::insert(hash_u hash)
{
    if ( use64 )
    {
        hashes64.insert(hash.hash64);
    }
    else
    {
        hashes32.insert(hash.hash32);
    }
}

void HashSet::toHashList(HashList & hashList) const
{
    if ( use64 )
    {
        for ( std::unordered_set<hash64_t>::const_iterator i = hashes64.begin(); i != hashes64.end(); i++ )
        {
            hashList.push_back64(*i);
        }
    }
    else
    {
        for ( std::unordered_set<hash32_t>::const_iterator i = hashes32.begin(); i != hashes32.end(); i++ )
        {
            hashList.push_back32(*i);
        }
    }
}
