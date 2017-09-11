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
#include <unordered_set>
#include <unordered_map>

namespace mash {

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
double estimateDistance(uint64_t common, uint64_t denom, int kmerSize, double kmerSpace);
double pValueWithin(uint64_t x, uint64_t setSize, double kmerSpace, uint64_t sketchSize);
void translate(const char * src, char * dst, uint64_t len);

} // namespace mash

#endif
