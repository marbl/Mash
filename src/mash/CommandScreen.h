// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandScreen
#define INCLUDED_CommandScreen

#include "Command.h"
#include "Sketch.h"
#include <list>
#include <string>
#include <vector>
#include <atomic>
#include <unordered_set>
#include <unordered_map>
#include "MinHashHeap.h"

namespace mash {

struct HashTableEntry
{
	HashTableEntry() : count(0) {}
	
	uint32_t count;
	std::unordered_set<uint64_t> indices;
};

//typedef std::unordered_map< uint64_t, HashTableEntry > HashTable;
typedef std::unordered_map< uint64_t, std::unordered_set<uint64_t> > HashTable;

static const std::unordered_map< std::string, char > codons =
{
	{"AAA",	'K'},
	{"AAC",	'N'},
	{"AAG",	'K'},
	{"AAT",	'N'},
	{"ACA",	'T'},
	{"ACC",	'T'},
	{"ACG",	'T'},
	{"ACT",	'T'},
	{"AGA",	'R'},
	{"AGC",	'S'},
	{"AGG",	'R'},
	{"AGT",	'S'},
	{"ATA",	'I'},
	{"ATC",	'I'},
	{"ATG",	'M'},
	{"ATT",	'I'},
	{"CAA",	'Q'},
	{"CAC",	'H'},
	{"CAG",	'Q'},
	{"CAT",	'H'},
	{"CCA",	'P'},
	{"CCC",	'P'},
	{"CCG",	'P'},
	{"CCT",	'P'},
	{"CGA",	'R'},
	{"CGC",	'R'},
	{"CGG",	'R'},
	{"CGT",	'R'},
	{"CTA",	'L'},
	{"CTC",	'L'},
	{"CTG",	'L'},
	{"CTT",	'L'},
	{"GAA",	'E'},
	{"GAC",	'D'},
	{"GAG",	'E'},
	{"GAT",	'D'},
	{"GCA",	'A'},
	{"GCC",	'A'},
	{"GCG",	'A'},
	{"GCT",	'A'},
	{"GGA",	'G'},
	{"GGC",	'G'},
	{"GGG",	'G'},
	{"GGT",	'G'},
	{"GTA",	'V'},
	{"GTC",	'V'},
	{"GTG",	'V'},
	{"GTT",	'V'},
	{"TAA",	'*'},
	{"TAC",	'Y'},
	{"TAG",	'*'},
	{"TAT",	'Y'},
	{"TCA",	'S'},
	{"TCC",	'S'},
	{"TCG",	'S'},
	{"TCT",	'S'},
	{"TGA",	'*'},
	{"TGC",	'C'},
	{"TGG",	'W'},
	{"TGT",	'C'},
	{"TTA",	'L'},
	{"TTC",	'F'},
	{"TTG",	'L'},
	{"TTT",	'F'}
};

class CommandScreen : public Command
{
public:
    
    struct HashInput
    {
    	HashInput(std::unordered_map<uint64_t, std::atomic<uint32_t> > & hashCountsNew, MinHashHeap * minHashHeapNew, char * seqNew, uint64_t lengthNew, const Sketch::Parameters & parametersNew, bool transNew)
    	:
    	hashCounts(hashCountsNew),
    	minHashHeap(minHashHeapNew),
    	seq(seqNew),
    	length(lengthNew),
    	parameters(parametersNew),
    	trans(transNew)
    	{}
    	
    	~HashInput()
    	{
    		if ( seq != 0 )
    		{
	    		delete [] seq;
	    	}
    	}
    	
    	std::string fileName;
    	
    	char * seq;
    	uint64_t length;
    	bool trans;
    	
    	Sketch::Parameters parameters;
		std::unordered_map<uint64_t, std::atomic<uint32_t> > & hashCounts;
		MinHashHeap * minHashHeap;
    };
    
    struct HashOutput
    {
    	HashOutput(MinHashHeap * minHashHeapNew)
    	:
    	minHashHeap(minHashHeapNew)
    	{}
    	
		MinHashHeap * minHashHeap;
    };
    
    CommandScreen();
    
    int run() const; // override

private:
	
	struct Reference
	{
		Reference(uint64_t amerCountNew, std::string nameNew, std::string commentNew)
		: amerCount(amerCountNew), name(nameNew), comment(commentNew) {}
		
		uint64_t amerCount;
		std::string name;
		std::string comment;
	};
};

char aaFromCodon(const char * codon);
double estimateIdentity(uint64_t common, uint64_t denom, int kmerSize, double kmerSpace);
CommandScreen::HashOutput * hashSequence(CommandScreen::HashInput * input);
double pValueWithin(uint64_t x, uint64_t setSize, double kmerSpace, uint64_t sketchSize);
void translate(const char * src, char * dst, uint64_t len);
void useThreadOutput(CommandScreen::HashOutput * output, std::unordered_set<MinHashHeap *> & minHashHeaps);

} // namespace mash

#endif
