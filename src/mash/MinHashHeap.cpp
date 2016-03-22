
#include "MinHashHeap.h"
#include <iostream>

using namespace::std;

MinHashHeap::MinHashHeap(bool use64New, uint64_t cardinalityMaximumNew, uint64_t multiplicityMinimumNew, uint64_t memoryBoundBytes) :
	use64(use64New),
	hashes(use64New),
	hashesQueue(use64New),
	hashesPending(use64New),
	hashesQueuePending(use64New)
{
	cardinalityMaximum = cardinalityMaximumNew;
	multiplicityMinimum = multiplicityMinimumNew;
	
	multiplicitySum = 0;
	
	if ( memoryBoundBytes == 0 )
	{
		bloomFilter = 0;
	}
	else
	{
		bloom_parameters bloomParams;
		
		bloomParams.projected_element_count = 1000000000;//(uint64_t)parameters.genomeSize * 10l; // TODO: error rate based on platform and coverage
		bloomParams.false_positive_probability = 0;//parameters.bloomError;
		bloomParams.maximum_size = memoryBoundBytes * 8l;
		bloomParams.compute_optimal_parameters();
		
		kmersTotal = 0;
		kmersUsed = 0;
		
		//if ( i == 0 && verbosity > 0 )
		{
			//cerr << "   Bloom table size (bytes): " << bloomParams.optimal_parameters.table_size / 8 << endl;
		}
		
		bloomFilter = new bloom_filter(bloomParams);
	}
}

MinHashHeap::~MinHashHeap()
{
	if ( bloomFilter != 0 )
	{
		delete bloomFilter;
	}
}

void MinHashHeap::computeStats()
{
	vector<uint32_t> counts;
	hashes.toCounts(counts);
	
	for ( int i = 0; i < counts.size(); i++ )
	{
		cout << counts.at(i) << endl;
	}
}

void MinHashHeap::clear()
{
	hashes.clear();
	hashesQueue.clear();
	
	hashesPending.clear();
	hashesQueuePending.clear();
	
	if ( bloomFilter != 0 )
	{
		bloomFilter->clear();
	}
	
	multiplicitySum = 0;
}

void MinHashHeap::tryInsert(hash_u hash)
{
	if
	(
		hashes.size() < cardinalityMaximum ||
		hashLessThan(hash, hashesQueue.top(), use64)
	)
	{
		if ( hashes.count(hash) == 0 )
		{
			if ( bloomFilter != 0 )
			{
                const unsigned char * data = use64 ? (const unsigned char *)&hash.hash64 : (const unsigned char *)&hash.hash32;
            	size_t length = use64 ? 8 : 4;
            	
                if ( bloomFilter->contains(data, length) )
                {
					hashes.insert(hash, 2);
					hashesQueue.push(hash);
					multiplicitySum += 2;
	                kmersUsed++;
                }
            	else
            	{
	                bloomFilter->insert(data, length);
	                kmersTotal++;
	            }
			}
			else if ( multiplicityMinimum == 1 || hashesPending.count(hash) == multiplicityMinimum - 1 )
			{
				hashes.insert(hash, multiplicityMinimum);
				hashesQueue.push(hash);
				multiplicitySum += multiplicityMinimum;
				
				if ( multiplicityMinimum > 1 )
				{
					// just remove from set for now; will be removed from
					// priority queue when it's on top
					//
					hashesPending.erase(hash);
				}
			}
			else
			{
				if ( hashesPending.count(hash) == 0 )
				{
					hashesQueuePending.push(hash);
				}
			
				hashesPending.insert(hash, 1);
			}
		}
		else
		{
			hashes.insert(hash, 1);
			multiplicitySum++;
		}
		
		if ( hashes.size() > cardinalityMaximum )
		{
			multiplicitySum -= hashes.count(hashesQueue.top());
			hashes.erase(hashesQueue.top());
			
			// loop since there could be zombie hashes (gone from hashesPending)
			//
			while ( hashesQueuePending.size() > 0 && hashLessThan(hashesQueue.top(), hashesQueuePending.top(), use64) )
			{
				if ( hashesPending.count(hashesQueuePending.top()) )
				{
					hashesPending.erase(hashesQueuePending.top());
				}
				
				hashesQueuePending.pop();
			}
			
			hashesQueue.pop();
		}
	}
}
