// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef HashList_h
#define HashList_h

#include "hash.h"
#include <vector>

class HashList
{
public:
    
    HashList() {use64 = true;}
    HashList(int kmerSize) {use64 = kmerSize > 16;}
    
    hash_u at(int index) const;
    void clear();
    void resize(int size);
    void set32(int index, uint32_t value);
    void set64(int index, uint64_t value);
    void setUse64(bool use64New) {use64 = use64New;}
    int size() const {return use64 ? hashes64.size() : hashes32.size();}
    void sort();
    void push_back32(hash32_t hash) {hashes32.push_back(hash);}
    void push_back64(hash64_t hash) {hashes64.push_back(hash);}
    bool get64() const {return use64;}
    
private:
    
    bool use64;
    std::vector<hash32_t> hashes32;
    std::vector<hash64_t> hashes64;
};

#endif
