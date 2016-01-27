#ifndef HashHeapCounted_h
#define HashHeapCounted_h

#include "HashList.h"
#include "HashPriorityQueue.h"
#include "HashSet.h"
#include <math.h>

class MinHashHeap
{
public:

	MinHashHeap(bool use64New, uint64_t cardinalityMaximumNew, uint64_t multiplicityMinimumNew = 1);
	void computeStats();
	void clear();
	double estimateMultiplicity() const;
	double estimateSetSize() const;
    void toHashList(HashList & hashList) const;
	void tryInsert(hash_u hash);

private:

	bool use64;
	
	HashSet hashes;
	HashPriorityQueue hashesQueue;
	
	HashSet hashesPending;
	HashPriorityQueue hashesQueuePending;
	
	uint64_t cardinalityMaximum;
	uint64_t multiplicityMinimum;
	
	uint64_t multiplicitySum;
};

inline double MinHashHeap::estimateMultiplicity() const {return (double)multiplicitySum / hashes.size();}
inline double MinHashHeap::estimateSetSize() const {return pow(2.0, use64 ? 64.0 : 32.0) * (double)hashes.size() / (use64 ? (double)hashesQueue.top().hash64 : (double)hashesQueue.top().hash32);}
inline void MinHashHeap::toHashList(HashList & hashList) const {hashes.toHashList(hashList);}

#endif
