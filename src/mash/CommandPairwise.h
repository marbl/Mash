// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandPairwise
#define INCLUDED_CommandPairwise

#include "Command.h"
#include "Sketch.h"
#include <list>
#include <vector>
#include <unordered_map>

//typedef std::list<uint64_t> * HashTable;

struct HashEntry
{
	HashEntry(uint64_t indexSeqNew, uint64_t indexHashNew) :
		indexSeq(indexSeqNew), indexHash(indexHashNew) {};
	
	uint64_t indexSeq;
	uint64_t indexHash;
};

struct Comparison
{
	Comparison()
	{
		shared = 0;
		skipped = 0;
		lastSharedIndexQry = 0;
		lastSharedIndexRef = 0;
	}
	
	uint64_t shared;
	uint64_t skipped;
	
	// 1-indexed for easier math
	//
	uint64_t lastSharedIndexQry;
	uint64_t lastSharedIndexRef;
};

typedef std::unordered_map< uint64_t, std::list<HashEntry> > HashTable;

class CommandPairwise : public Command
{
public:
    
    struct PairwiseInput
    {
        PairwiseInput(const Sketch & sketchNew, uint64_t indexNew, const Sketch::Parameters & parametersNew, double maxDistanceNew, double maxPValueNew, const HashTable & hashTableNew, uint64_t hashTableSizeNew)
            :
            sketch(sketchNew),
            index(indexNew),
            parameters(parametersNew),
            maxDistance(maxDistanceNew),
            maxPValue(maxPValueNew),
            hashTable(hashTableNew),
            hashTableSize(hashTableSizeNew)
            {}
        
        const Sketch & sketch;
        
        uint64_t index;
        
        const Sketch::Parameters & parameters;
        double maxDistance;
        double maxPValue;
        
		const HashTable & hashTable;
		uint64_t hashTableSize;
    };
    
    struct PairwiseOutput
    {
        PairwiseOutput(const Sketch & sketchNew, uint64_t indexNew)
            :
            sketch(sketchNew),
            index(indexNew)
        {
        }
        
        ~PairwiseOutput()
        {
        }
        
        struct PairOutput
        {
        	uint64_t index;
            uint64_t numer;
            uint64_t denom;
            double distance;
            double pValue;
        };
        
        static bool pairOutputLessThan(const PairOutput & a, const PairOutput & b);
        
        const Sketch & sketch;
        
        uint64_t index;
        
        std::vector<PairOutput> pairs;
    };
    
    CommandPairwise();
    
    int run() const; // override
    
private:
    
    void writeOutput(PairwiseOutput * output, bool table) const;
};

void fillHashTable(const Sketch & sketch, HashTable & hashTable, uint64_t start, uint64_t end);
CommandPairwise::PairwiseOutput * search(CommandPairwise::PairwiseInput * input);
bool compareSketches(CommandPairwise::PairwiseOutput::PairOutput * output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t common, uint64_t denom, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxDistance, double maxPValue);

#endif
