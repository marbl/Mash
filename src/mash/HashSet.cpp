// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "HashSet.h"
#include <utility>

uint32_t HashSet::count(hash_u hash) const
{
	if ( use64 )
	{
		if ( hashes64.count(hash.hash64) )
		{
			return hashes64.at(hash.hash64);
		}
		else
		{
			return 0;
		}
	}
	else
	{
		if ( hashes32.count(hash.hash32) )
		{
			return hashes32.at(hash.hash32);
		}
		else
		{
			return 0;
		}
	}
}

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

void HashSet::insert(hash_u hash, uint32_t count)
{
    if ( use64 )
    {
    	hash64_t hash64 = hash.hash64;
    	
    	if ( hashes64.count(hash64) )
    	{
    		hashes64[hash64] = hashes64.at(hash64) + count;
    	}
    	else
    	{
	        hashes64[hash64] = count;
	    }
    }
    else
    {
    	hash32_t hash32 = hash.hash32;
    	
    	if ( hashes32.count(hash32) )
    	{
    		hashes32[hash32] = hashes32.at(hash32) + count;
    	}
    	else
    	{
	        hashes32[hash32] = count;
	    }
    }
}

void HashSet::toHashList(HashList & hashList, std::vector<uint32_t> & counts) const
{
    if ( use64 )
    {
    	typedef std::pair<hash64_t, uint32_t> hashCount64;
    	
    	std::vector <hashCount64> sortList;
    	
        for ( auto i = hashes64.begin(); i != hashes64.end(); i++ )
        {
            sortList.push_back(hashCount64(i->first, i->second));
        }
        
        std::sort(sortList.begin(), sortList.end(), [](hashCount64 a, hashCount64 b){return a.first < b.first;});
        
        for ( auto i = sortList.begin(); i != sortList.end(); i++ )
        {
            hashList.push_back64(i->first);
            counts.push_back(i->second);
        }
    }
    else
    {
    	typedef std::pair<hash32_t, uint32_t> hashCount32;
    	
    	std::vector <hashCount32> sortList;
    	
        for ( auto i = hashes32.begin(); i != hashes32.end(); i++ )
        {
            sortList.push_back(hashCount32(i->first, i->second));
        }
        
        std::sort(sortList.begin(), sortList.end(), [](hashCount32 a, hashCount32 b){return a.first < b.first;});
        
        for ( auto i = sortList.begin(); i != sortList.end(); i++ )
        {
            hashList.push_back32(i->first);
            counts.push_back(i->second);
        }
    }
}

void HashSet::toHashList(HashList & hashList) const
{
    if ( use64 )
    {
        for ( auto i = hashes64.begin(); i != hashes64.end(); i++ )
        {
            hashList.push_back64(i->first);
        }
    }
    else
    {
        for ( auto i = hashes32.begin(); i != hashes32.end(); i++ )
        {
            hashList.push_back32(i->first);
        }
    }
}
