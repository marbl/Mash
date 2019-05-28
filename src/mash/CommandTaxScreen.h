// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandTaxScreen
#define INCLUDED_CommandTaxScreen

#include "Command.h"
#include "Sketch.h"
#include <list>
#include <string>
#include <vector>
#include <atomic>
#include <unordered_set>
#include <unordered_map>
#include "MinHashHeap.h"
#include "CommandScreen.h"


using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::list;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;


namespace mash {

using TaxID = uint64_t;

typedef std::unordered_map< uint64_t, std::unordered_set<uint64_t> > HashTable;

class CommandTaxScreen : public Command
{
public:
    
    CommandTaxScreen();
    
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

} // namespace mash

#endif
