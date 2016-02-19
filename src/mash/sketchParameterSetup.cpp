// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "sketchParameterSetup.h"
#include <iostream>

using namespace::std;

int sketchParameterSetup(Sketch::Parameters & parameters, const Command & command)
{
    parameters.kmerSize = command.getOption("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = command.getOption("sketchSize").getArgumentAsNumber();
    parameters.concatenated = ! command.getOption("individual").active;
    parameters.noncanonical = command.getOption("noncanonical").active;
    parameters.reads = command.getOption("reads").active;
    parameters.minCov = command.getOption("minCov").getArgumentAsNumber();
    parameters.targetCov = command.getOption("targetCov").getArgumentAsNumber();
    parameters.windowed = false;//command.getOption("windowed").active;
    parameters.windowSize = 0;//command.getOption("window").getArgumentAsNumber();
    parameters.parallelism = command.getOption("threads").getArgumentAsNumber();
    parameters.preserveCase = command.getOption("case").active;
    
	if ( command.hasOption("warning") )
	{
	    parameters.warning = command.getOption("warning").getArgumentAsNumber();
	}
	
    if ( command.getOption("minCov").active || command.getOption("targetCov").active )
    {
        parameters.reads = true;
    }
    
    if ( parameters.reads && ! parameters.concatenated )
    {
        cerr << "ERROR: The option " << command.getOption("individual").identifier << " cannot be used with " << command.getOption("unique").identifier << "." << endl;
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

