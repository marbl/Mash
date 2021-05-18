// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "sketchParameterSetup.h"
#include <iostream>

using std::cerr;
using std::endl;

namespace mash {

int sketchParameterSetup(Sketch::Parameters & parameters, const Command & command)
{
    parameters.kmerSize = command.getOption("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = command.getOption("sketchSize").getArgumentAsNumber();
    parameters.concatenated = ! command.getOption("individual").active;
    parameters.noncanonical = command.getOption("noncanonical").active;
    parameters.seed = command.getOption("seed").getArgumentAsNumber();
    parameters.reads = command.getOption("reads").active;
    parameters.minCov = command.getOption("minCov").getArgumentAsNumber();
    parameters.targetCov = command.getOption("targetCov").getArgumentAsNumber();
#ifdef COMMAND_FIND
    parameters.windowed = command.getOption("windowed").active;
    parameters.windowSize = command.getOption("window").getArgumentAsNumber();
    parameters.concatenated = false;
#endif
    parameters.parallelism = command.getOption("threads").getArgumentAsNumber();
    parameters.preserveCase = command.getOption("case").active;
    
	if ( command.hasOption("warning") )
	{
	    parameters.warning = command.getOption("warning").getArgumentAsNumber();
	}
	
	if ( command.getOption("memory").active )
	{
		parameters.reads = true;
		//parameters.minCov = 2;
		parameters.memoryBound = command.getOption("memory").getArgumentAsNumber();
		
		if ( command.getOption("minCov").active )
		{
			cerr << "ERROR: The option " << command.getOption("minCov").identifier << " cannot be used with " << command.getOption("memory").identifier << "." << endl;
			return 1;
		}
	}
	
    if ( command.getOption("minCov").active || command.getOption("targetCov").active )
    {
        parameters.reads = true;
    }
    
    if ( command.getOption("genome").active )
    {
    	parameters.reads = true;
    	parameters.genomeSize = command.getOption("genome").getArgumentAsNumber();
    }
    
    if ( parameters.reads )
    {
    	parameters.counts = true;
    }
    
    if ( parameters.reads && command.getOption("threads").active )
    {
    	cerr << "WARNING: The option " << command.getOption("threads").identifier << " will be ignored with " << command.getOption("reads").identifier << "." << endl;
    }
    
    if ( parameters.reads && ! parameters.concatenated )
    {
        cerr << "ERROR: The option " << command.getOption("individual").identifier << " cannot be used with " << command.getOption("reads").identifier << "." << endl;
        return 1;
    }
    
    if ( parameters.concatenated && parameters.windowed )
    {
        cerr << "ERROR: " << command.getOption("concat").identifier << " and " << command.getOption("windowed").identifier << " are incompatible." << endl;
        return 1;
    }
    
    if ( command.getOption("protein").active )
    {
    	parameters.noncanonical = true;
    	setAlphabetFromString(parameters, alphabetProtein);
    	
    	if ( ! command.getOption("kmer").active )
    	{
    		parameters.kmerSize = 9;
    	}
    }
    else if ( command.getOption("alphabet").active )
    {
    	parameters.noncanonical = true;
    	setAlphabetFromString(parameters, command.getOption("alphabet").argument.c_str());
    }
    else
    {
    	setAlphabetFromString(parameters, alphabetNucleotide);
    }
    
    return 0;
}

void warnKmerSize(const Sketch::Parameters & parameters, const Command & command, uint64_t lengthMax, const std::string & lengthMaxName, double randomChance, int kMin, int warningCount)
{
	cerr << "\nWARNING: For the k-mer size used (" << parameters.kmerSize
		<< "), the random match probability (" << randomChance
		<< ") is above the specified warning threshold ("
		<< parameters.warning << ") for the sequence \"" << lengthMaxName
		<< "\" of size " << lengthMax;
	
	if ( warningCount > 1 )
	{
		cerr << " (and " << (warningCount - 1) << " others)";
	}
	
	cerr << ". Distances to "
		<< (warningCount == 1 ? "this sequence" : "these sequences")
		<< " may be underestimated as a result. To meet the threshold of "
		<< parameters.warning << ", a k-mer size of at least " << kMin
		<< " is required. See: -" << command.getOption("kmer").identifier << ", -" << command.getOption("warning").identifier << "." << endl << endl;
}

} // namespace mash