// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandPairwise.h"
#include "CommandDistance.h" // for pvalue
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"
#include "sketchParameterSetup.h"
#include <math.h>
#include <set>

#ifdef USE_BOOST
	#include <boost/math/distributions/binomial.hpp>
	using namespace::boost::math;
#else
	#include <gsl/gsl_cdf.h>
#endif

using namespace::std;

bool CommandPairwise::PairwiseOutput::pairOutputLessThan(const PairOutput & a, const PairOutput & b)
{
	return a.index < b.index;
}

CommandPairwise::CommandPairwise()
: Command()
{
	name = "pairwise";
	summary = "Estimate the pairwise distance of protein sequences.";
	description = "Estimate the pairwise distance of protein sequences. Input can be fasta or a Mash sketch file (.msh). The output fields are [seq-ID-1, seq-ID-2, distance, p-value, shared-hashes].";
	argumentString = "<fasta>";
	
	useOption("help");
	useOption("threads");
    addOption("kmer", Option(Option::Integer, "k", "Sketch", "K-mer size. Hashes will be based on strings of this many amino acids.", "7", 1, 32));
    addOption("sketchSize", Option(Option::Integer, "s", "Sketch", "Sketch size. Each sketch will have at most this many non-redundant min-hashes.", "400"));
    addOption("case", Option(Option::Boolean, "Z", "Sketch", "Preserve case in k-mers and alphabet (case is ignored by default). Sequence letters whose case is not in the current alphabet will be skipped when sketching.", ""));
//	addOption("list", Option(Option::Boolean, "l", "Input", "List input. Each query file contains a list of sequence files, one per line. The reference file is not affected.", ""));
	addOption("table", Option(Option::Boolean, "t", "Output", "Table output (will not report p-values, but fields will be blank if they do not meet the p-value threshold).", ""));
	//addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances and divide by k-mer size to provide a better analog to phylogenetic distance. The special case of zero shared min-hashes will result in a distance of 1.", ""));
	addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
	addOption("distance", Option(Option::Number, "d", "Output", "Maximum distance to report.", "1.0", 0., 1.));
	//useSketchOptions();
}

int CommandPairwise::run() const
{
	if ( arguments.size() != 1 || options.at("help").active )
	{
		print();
		return 0;
	}
	
	int threads = options.at("threads").getArgumentAsNumber();
	bool table = options.at("table").active;
	//bool log = options.at("log").active;
	double pValueMax = options.at("pvalue").getArgumentAsNumber();
	double distanceMax = options.at("distance").getArgumentAsNumber();
	
	Sketch::Parameters parameters;
	
    parameters.kmerSize = getOption("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = getOption("sketchSize").getArgumentAsNumber();
    parameters.parallelism = getOption("threads").getArgumentAsNumber();
    parameters.preserveCase = getOption("case").active;
	parameters.noncanonical = true;
	parameters.concatenated = false;
	setAlphabetFromString(parameters, alphabetProtein);
	Sketch sketchRef;
	
	uint64_t lengthMax;
	double randomChance;
	int kMin;
	string lengthMaxName;
	int warningCount = 0;
	
	const string & fileReference = arguments[0];
	
	bool isSketch = hasSuffix(fileReference, suffixSketch);
	
	if ( isSketch )
	{
		if ( options.at("kmer").active )
		{
			cerr << "ERROR: The option -" << options.at("kmer").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
			return 1;
		}
		
		if ( options.at("noncanonical").active )
		{
			cerr << "ERROR: The option -" << options.at("noncanonical").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
			return 1;
		}
		
		if ( options.at("protein").active )
		{
			cerr << "ERROR: The option -" << options.at("protein").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
			return 1;
		}
		
		if ( options.at("alphabet").active )
		{
			cerr << "ERROR: The option -" << options.at("alphabet").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
			return 1;
		}
	}
	else
	{
		cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...";
	}
	
	vector<string> refArgVector;
	refArgVector.push_back(fileReference);
	
	//cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;
	
	sketchRef.initFromFiles(refArgVector, parameters);
	
	double lengthThreshold = (parameters.warning * sketchRef.getKmerSpace()) / (1. - parameters.warning);
	
	if ( isSketch )
	{
		if ( options.at("sketchSize").active )
		{
			if ( parameters.reads && parameters.minHashesPerWindow != sketchRef.getMinHashesPerWindow() )
			{
				cerr << "ERROR: The sketch size must match the reference when using a bloom filter (leave this option out to inherit from the reference sketch)." << endl;
				return 1;
			}
		}
		
		parameters.minHashesPerWindow = sketchRef.getMinHashesPerWindow();
		parameters.kmerSize = sketchRef.getKmerSize();
		parameters.noncanonical = sketchRef.getNoncanonical();
		parameters.preserveCase = sketchRef.getPreserveCase();
		
		string alphabet;
		sketchRef.getAlphabetAsString(alphabet);
		setAlphabetFromString(parameters, alphabet.c_str());
	}
	else
	{
		for ( uint64_t i = 0; i < sketchRef.getReferenceCount(); i++ )
		{
			uint64_t length = sketchRef.getReference(i).length;
		
			if ( length > lengthThreshold )
			{
				if ( warningCount == 0 || length > lengthMax )
				{
					lengthMax = length;
					lengthMaxName = sketchRef.getReference(i).name;
					randomChance = sketchRef.getRandomKmerChance(i);
					kMin = sketchRef.getMinKmerSize(i);
				}
			
				warningCount++;
			}
		}
	
		cerr << "done.\n";
	}
	
	if ( table )
	{
		cout << "#query";
		
		for ( int i = 0; i < sketchRef.getReferenceCount(); i++ )
		{
			cout << '\t' << sketchRef.getReference(i).name;
		}
		
		cout << endl;
	}
	
	ThreadPool<PairwiseInput, PairwiseOutput> threadPool(search, threads);
	
	vector<string> queryFiles;
	
	for ( int i = 0; i < arguments.size(); i++ )
	{
		queryFiles.push_back(arguments[i]);
	}
	
	Sketch sketch;
	
	sketch.initFromFiles(queryFiles, parameters);
	
	HashTable hashTable;
	uint64_t hashTableSize = fillHashTable(sketch, hashTable);
	
	for ( uint64_t i = 1; i < sketch.getReferenceCount(); i++ )
	{
		threadPool.runWhenThreadAvailable(new PairwiseInput(sketch, i, parameters, distanceMax, pValueMax, hashTable, hashTableSize));
		
		while ( threadPool.outputAvailable() )
		{
			writeOutput(threadPool.popOutputWhenAvailable(), table);
		}
	}
	
	while ( threadPool.running() )
	{
		writeOutput(threadPool.popOutputWhenAvailable(), table);
	}
	
	if ( warningCount > 0 && ! parameters.reads )
	{
		//warnKmerSize(parameters, *this, lengthMax, lengthMaxName, randomChance, kMin, warningCount);
	}
	
	//delete [] hashTable;
	
	return 0;
}

void CommandPairwise::writeOutput(PairwiseOutput * output, bool table) const
{
	uint64_t i = output->index;
	uint64_t j = 0;
	
	if ( table )
	{
		cout << output->sketch.getReference(output->index).name << '\t';
	}
	
	for ( uint64_t p = 0; p < output->pairs.size(); p++ )
	{
		const PairwiseOutput::PairOutput * pair = &output->pairs.at(p);
		
		//cerr << "index: " << pair->index << endl;
		if ( table )
		{
			while ( j < pair->index )
			{
				cout << '\t';
				j++;
			}
			
			cout << '\t' << pair->distance;
		}
		else
		{
			cout << output->sketch.getReference(i).name << '\t' << output->sketch.getReference(pair->index).name << '\t' << pair->distance << '\t' << pair->pValue << '\t' << pair->numer << '/' << pair->denom << endl;
		}
	}
	
	if ( table )
	{
		cout << endl;
	}
	
	delete output;
}

uint64_t fillHashTable(const Sketch & sketch, HashTable & hashTable)
{
	uint64_t hashTableSize = sketch.getKmerSpace();
	//hashTable = new list<uint64_t>[hashTableSize];
	
	cerr << "Creating hash table..." << endl;
	for ( int i = 0; i < sketch.getReferenceCount(); i++ )
	{
//		cerr << "i: " << i << endl;
		const HashList & hashesSorted = sketch.getReference(i).hashesSorted;
		
		for ( uint64_t j = 0; j < hashesSorted.size(); j++ )
		{
//			cerr << "hash: " << 1 << endl;
			hashTable[hashesSorted.at(j).hash32].push_back(i);
		}
	}
	
	return hashTableSize;
}

CommandPairwise::PairwiseOutput * search(CommandPairwise::PairwiseInput * input)
{
	const Sketch & sketch = input->sketch;
	
	CommandPairwise::PairwiseOutput * output = new CommandPairwise::PairwiseOutput(input->sketch, input->index);
	
	uint64_t sketchSize = sketch.getMinHashesPerWindow();
	
	const HashTable & hashTable = input->hashTable;
	uint64_t hashTableSize = input->hashTableSize;
	
	const HashList & hashesSortedRef = sketch.getReference(input->index).hashesSorted;
	
	set<uint64_t> targets;
	
	for ( uint64_t i = 0; i < hashesSortedRef.size(); i++ )
	{
		//cout << "hash: " << (hashesSortedRef.at(i).hash32 % hashTableSize) << endl;
		const list<uint64_t> & indeces = hashTable.at(hashesSortedRef.at(i).hash32);
		
		for ( list<uint64_t>::const_iterator j = indeces.begin(); j != indeces.end(); j++ )
		{
			uint64_t index = *j;
			//cout << "  index: " << index << endl;
			if ( index < input->index )
			{
				targets.insert(index);
			}
		}
	}
	
	for ( set<uint64_t>::const_iterator i = targets.begin(); i != targets.end(); i++ )
	{
		CommandPairwise::PairwiseOutput::PairOutput pair;
		
//		cerr << "comparing " << input->index << " to " << *i << endl;
		if ( compareSketches(&pair, sketch.getReference(input->index), sketch.getReference(*i), sketchSize, sketch.getKmerSize(), sketch.getKmerSpace(), input->maxDistance, input->maxPValue) )
		{
			//cerr << "hit!" << endl;
			pair.index = *i;
			output->pairs.push_back(pair);
		}
	}
	
	//sort(output->pairs.begin(), output->pairs.end(), CommandPairwise::PairwiseOutput::pairOutputLessThan);
	
	return output;
}

bool compareSketches(CommandPairwise::PairwiseOutput::PairOutput * output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxDistance, double maxPValue)
{
	uint64_t i = 0;
	uint64_t j = 0;
	uint64_t common = 0;
	uint64_t denom = 0;
	const HashList & hashesSortedRef = refRef.hashesSorted;
	const HashList & hashesSortedQry = refQry.hashesSorted;
	
	while ( denom < sketchSize && i < hashesSortedRef.size() && j < hashesSortedQry.size() )
	{
		//cerr << hashesSortedRef.at(i).hash32 << '\t' << hashesSortedQry.at(j).hash32 << endl;
		if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQry.at(j), hashesSortedRef.get64()) )
		{
			i++;
		}
		else if ( hashLessThan(hashesSortedQry.at(j), hashesSortedRef.at(i), hashesSortedRef.get64()) )
		{
			j++;
		}
		else
		{
			i++;
			j++;
			common++;
		}
		
		denom++;
	}
	
	if ( denom < sketchSize )
	{
		// complete the union operation if possible
		
		if ( i < hashesSortedRef.size() )
		{
			denom += hashesSortedRef.size() - i;
		}
		
		if ( j < hashesSortedQry.size() )
		{
			denom += hashesSortedQry.size() - j;
		}
		
		if ( denom > sketchSize )
		{
			denom = sketchSize;
		}
	}
	
	double distance;
	double jaccard = double(common) / denom;
	
	if ( common == denom ) // avoid -0
	{
		distance = 0;
	}
	else if ( common == 0 ) // avoid inf
	{
		distance = 1.;
	}
	else
	{
		//distance = log(double(common + 1) / (denom + 1)) / log(1. / (denom + 1));
		distance = -log(2 * jaccard / (1. + jaccard)) / kmerSize;
	}
	
	if ( distance > maxDistance || distance == 1 )
	{
		return false;
	}
	
	output->numer = common;
	output->denom = denom;
	output->distance = distance;
	output->pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);
	//cerr << "sketch size: " << sketchSize << " common: " << common << " denom: " << denom << " dist: " << distance << " pval: " << output->pValue <<  endl;
	
	return output->pValue <= maxPValue;
}

