// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef HashSet_h
#define HashSet_h

#include "HashList.h"
#include <unordered_set>

class HashSet
{
public:

    HashSet(int kmerSize) {use64 = kmerSize > 16;}
    
    int size() const {return use64 ? hashes64.size() : hashes32.size();}
    void clear() {use64 ? hashes64.clear() : hashes32.clear();}
    bool contains(hash_u hash) const {return use64 ? hashes64.count(hash.hash64) : hashes32.count(hash.hash32);}
    void erase(hash_u hash);
    void insert(hash_u hash);
    void toHashList(HashList & hashList) const;
    
private:
    
    bool use64;
    std::unordered_set<hash32_t> hashes32;
    std::unordered_set<hash64_t> hashes64;
};

#endif
