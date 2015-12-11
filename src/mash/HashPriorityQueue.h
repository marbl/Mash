// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef HashPriorityQueue_h
#define HashPriorityQueue_h

#include "hash.h"
#include <queue>

class HashPriorityQueue
{
public:
	
	HashPriorityQueue(bool use64New) : use64(use64New) {}
	void clear();
	void pop() {use64 ? queue64.pop() : queue32.pop();}
	void push(hash_u hash) {use64 ? queue64.push(hash.hash64) : queue32.push(hash.hash32);}
	int size() const {return use64 ? queue64.size() : queue32.size();}
	hash_u top() const;
	
private:
    
	bool use64;
	std::priority_queue<hash32_t> queue32;
	std::priority_queue<hash64_t> queue64;
};

#endif
