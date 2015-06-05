#include "HashList.h"
#include <algorithm>

hash_u HashList::at(int index) const
{
    hash_u hash;
    
    if ( use64 )
    {
        hash.hash64 = hashes64.at(index);
    }
    else
    {
        hash.hash32 = hashes32.at(index);
    }
    
    return hash;
}

void HashList::clear()
{
    if ( use64 )
    {
        hashes64.clear();
    }
    else
    {
        hashes32.clear();
    }
}

void HashList::resize(int size)
{
    if ( use64 )
    {
        hashes64.resize(size);
    }
    else
    {
        hashes32.resize(size);
    }
}

void HashList::set32(int index, uint32_t value)
{
    hashes32[index] = value;
}

void HashList::set64(int index, uint64_t value)
{
    hashes64[index] = value;
}

void HashList::sort()
{
    if ( use64 )
    {
        std::sort(hashes64.begin(), hashes64.end());
    }
    else
    {
        std::sort(hashes32.begin(), hashes32.end());
    }
}
