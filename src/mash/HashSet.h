// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef HashSet_h
#define HashSet_h

#include "HashList.h"
#include "robin_hood.h"
#include <vector>

class HashSet
{
public:

    HashSet(bool use64New) : use64(use64New) {}
    
    int size() const {return use64 ? hashes64.size() : hashes32.size();}
    void clear() {use64 ? hashes64.clear() : hashes32.clear();}
    uint32_t count(hash_u hash) const;
    void erase(hash_u hash);
    void insert(hash_u hash, uint32_t count = 1);
    void toHashList(HashList & hashList) const;
    void toCounts(std::vector<uint32_t> & counts) const;
    
private:
    
    bool use64;
    robin_hood::unordered_map<hash32_t, uint32_t> hashes32;
    robin_hood::unordered_map<hash64_t, uint32_t> hashes64;
};

#endif
