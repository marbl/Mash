// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandGenes.h"
#include "CommandDistance.h" // for pvalue
#include "Sketch.h"
#include "kseq.h"
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

#define SET_BINARY_MODE(file)
KSEQ_INIT(gzFile, gzread)

using namespace::std;

CommandGenes::CommandGenes()
: Command()
{
	name = "x";
	summary = "Estimate the pairwise distance of protein sequences.";
	description = "Estimate the pairwise distance of protein sequences. Input can be fasta or a Mash sketch file (.msh). The output fields are [seq-ID-1, seq-ID-2, distance, p-value, shared-hashes].";
	argumentString = "<fasta>";
	
	useOption("help");
	useOption("threads");
	useOption("minCov");
	useOption("targetCov");
    addOption("kmer", Option(Option::Integer, "k", "Sketch", "K-mer size. Hashes will be based on strings of this many amino acids.", "9", 1, 14));
    addOption("sketchSize", Option(Option::Integer, "s", "Sketch", "Sketch size. Each sketch will have at most this many non-redundant min-hashes.", "400"));
    addOption("case", Option(Option::Boolean, "Z", "Sketch", "Preserve case in k-mers and alphabet (case is ignored by default). Sequence letters whose case is not in the current alphabet will be skipped when sketching.", ""));
//	addOption("list", Option(Option::Boolean, "l", "Input", "List input. Each query file contains a list of sequence files, one per line. The reference file is not affected.", ""));
	addOption("table", Option(Option::Boolean, "t", "Output", "Table output (will not report p-values, but fields will be blank if they do not meet the p-value threshold).", ""));
	//addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances and divide by k-mer size to provide a better analog to phylogenetic distance. The special case of zero shared min-hashes will result in a distance of 1.", ""));
	addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
	addOption("distance", Option(Option::Number, "d", "Output", "Maximum distance to report.", "1.0", 0., 1.));
	//useSketchOptions();
}

int CommandGenes::run() const
{
	if ( arguments.size() < 2 || options.at("help").active )
	{
		print();
		return 0;
	}
	
	int threads = options.at("threads").getArgumentAsNumber();
	bool table = options.at("table").active;
	//bool log = options.at("log").active;
	double pValueMax = options.at("pvalue").getArgumentAsNumber();
	double distanceMax = options.at("distance").getArgumentAsNumber();
	int amerSize = getOption("kmer").getArgumentAsNumber();
	
	// for estimating saturation in nucleotide space
	//
	Sketch::Parameters parameters;
	//
	parameters.kmerSize = amerSize * 3;
	parameters.minHashesPerWindow = getOption("sketchSize").getArgumentAsNumber();
	parameters.parallelism = 1;
	parameters.preserveCase = getOption("case").active;
	parameters.noncanonical = false;
	parameters.concatenated = false;
    parameters.minCov = getOption("minCov").getArgumentAsNumber();
    parameters.targetCov = getOption("targetCov").getArgumentAsNumber();
	setAlphabetFromString(parameters, alphabetNucleotide);
	
	uint64_t amerSpace = pow(20, amerSize);
    parameters.use64 = amerSpace > pow(2, 32);
	
	const string & fileReference = arguments[0];
	
	//cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;
	
	HashTable hashTable;
	unordered_map<string, uint32_t> amerCounts;
	vector<CommandGenes::Reference> references;
	
	gzFile fp = gzopen(arguments[0].c_str(), "r");
	kseq_t * kseq = kseq_init(fp);
	
	cerr << "Filling table from " << arguments[0] << endl;
	
	int l;
	
	while ((l = kseq_read(kseq)) >= 0)
	{
		if ( l < amerSize )
		{
			continue;
		}
		
		for ( int i = 0; i < l - amerSize + 1; i++ )
		{
			string amer(kseq->seq.s + i, amerSize);
			hashTable[amer].insert(references.size());
			amerCounts[amer] = 0;
		}
		
		references.push_back(Reference(l - amerSize + 1, kseq->name.s, kseq->comment.s));
	}
	
	kseq_destroy(kseq);
	gzclose(fp);
	
	uint64_t * shared = new uint64_t[references.size()];
	
	memset(shared, 0, sizeof(uint64_t) * references.size());
	
	MinHashHeap minHashHeap(parameters.use64, parameters.minHashesPerWindow, parameters.reads ? parameters.minCov : 1, parameters.memoryBound);
	
	int queryCount = arguments.size() - 1;
	cerr << "Translating from " << queryCount << " inputs..." << endl;
	
	// open all query files for round robin
	//
	gzFile fps[queryCount];
	list<kseq_t *> kseqs;
	//
	for ( int f = 1; f < arguments.size(); f++ )
	{
		if ( arguments[f] == "-" )
		{
			if ( f > 1 )
			{
				cerr << "ERROR: '-' for stdin must be first query" << endl;
				exit(1);
			}
			
			fps[f - 1] = gzdopen(fileno(stdin), "r");
		}
		else
		{
			fps[f - 1] = gzopen(arguments[f].c_str(), "r");
		}
		
		kseqs.push_back(kseq_init(fps[f - 1]));
	}
	
	// perform round-robin, closing files as they end
	//
	uint64_t count = 0;
	list<kseq_t *>::iterator it = kseqs.begin();
	//
	while ( kseqs.begin() != kseqs.end() )
	{
		l = kseq_read(*it);
		
		if ( l < -1 ) // error
		{
			break;
		}
		
		if ( l == -1 ) // eof
		{
			kseq_destroy(*it);
			it = kseqs.erase(it);
			continue;
		}
		
		if ( l < amerSize ) // too short
		{
			continue;
		}
		
		count++;
		
		char * seq = (*it)->seq.s;
		
		it++;
		
		if ( it == kseqs.end() )
		{
			it = kseqs.begin();
		}
		
		// uppercase
		//
		for ( uint64_t i = 0; i < l; i++ )
		{
			if ( ! parameters.preserveCase && seq[i] > 96 && seq[i] < 123 )
			{
				seq[i] -= 32;
			}
		}
		
		char * seqRev = new char[l];
		
		reverseComplement(seq, seqRev, l);
		
		for ( int i = 0; i < 6; i++ )
		{
			int frame = i % 3;
			bool rev = i > 2;
			
			int lenTrans = (l - frame) / 3;
			
			char * seqTrans = new char[lenTrans];
			
			translate(rev ? seqRev : seq, seqTrans, lenTrans);
			
			string strTrans(seqTrans, lenTrans);
			//cout << strTrans << endl;
			
			int64_t lastGood = -1;
			
			for ( int j = 0; j < lenTrans - amerSize + 1; j++ )
			{
				while ( lastGood < j + amerSize && lastGood < lenTrans )
				{
					lastGood++;
					
					if ( seqTrans[lastGood + 1] == 0 )
					{
						j = lastGood + 1;
					}
				}
				
				if ( j > lenTrans - amerSize )
				{
					break;
				}
				
				string amer(seqTrans + j, amerSize);
				
				//cout << amer << endl;
				
				if ( hashTable.count(amer) == 1 )
				{
					amerCounts[amer]++;
					
					if ( amerCounts.at(amer) == parameters.minCov )
					{
						const unordered_set<uint64_t> & indeces = hashTable.at(amer);
				
						for ( unordered_set<uint64_t>::const_iterator k = indeces.begin(); k != indeces.end(); k++ )
						{
							shared[*k]++;
						}
					}
				}
			}
			
			delete [] seqTrans;
		}
		
		delete [] seqRev;
		
		addMinHashes(minHashHeap, seq, l, parameters);
		
		if ( parameters.targetCov > 0 && minHashHeap.estimateMultiplicity() >= parameters.targetCov )
		{
			l = -1; // success code
			break;
		}
	}
	
	for ( int i = 0; i < queryCount; i++ )
	{
		gzclose(fps[i]);
	}
	
	if (  l != -1 )
	{
		cerr << "\nERROR: reading inputs" << endl;
		exit(1);
	}
	
	if ( count == 0 )
	{
		cerr << "\nERROR: Did not find sequence records in inputs" << endl;
		
		exit(1);
	}
	
	cerr << "Reads required for " << parameters.targetCov << "x coverage: " << count << endl;
	cerr << "Estimated genome size: " << minHashHeap.estimateSetSize() << endl;
	
	for ( int i = 0; i < references.size(); i++ )
	{
		if ( shared[i] != 0 )
		{
			double identity = estimateIdentity(shared[i], references[i].amerCount, amerSize, amerSpace);
			
			cout << references[i].name << '\t' << identity << endl;
		}
	}
	
	delete [] shared;
	
	return 0;
}

double estimateIdentity(uint64_t common, uint64_t denom, int kmerSize, double kmerSpace)
{
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
	
	return 1. - distance;
}

void translate(const char * src, char * dst, uint64_t len)
{
	for ( uint64_t n = 0, a = 0; a < len; a++, n+= 3 )
	{
		dst[a] = aaFromCodon(src + n);
	}
}

char aaFromCodon(const char * codon)
{
	string str(codon, 3);
	
	if ( codons.count(str) == 1 )
	{
		return codons.at(str);
	}
	else
	{
		return 0;
	}
	char aa = 0;
	
	switch (codon[0])
	{
	case 'A':
		switch (codon[1])
		{
		case 'A':
			switch (codon[2])
			{
				case 'A': aa = 'K'; break;
				case 'C': aa = 'N'; break;
				case 'G': aa = 'K'; break;
				case 'T': aa = 'N'; break;
			}
			break;
		case 'C':
			switch (codon[2])
			{
				case 'A': aa = 'T'; break;
				case 'C': aa = 'T'; break;
				case 'G': aa = 'T'; break;
				case 'T': aa = 'T'; break;
			}
			break;
		case 'G':
			switch (codon[2])
			{
				case 'A': aa = 'R'; break;
				case 'C': aa = 'S'; break;
				case 'G': aa = 'R'; break;
				case 'T': aa = 'S'; break;
			}
			break;
		case 'T':
			switch (codon[2])
			{
				case 'A': aa = 'I'; break;
				case 'C': aa = 'I'; break;
				case 'G': aa = 'M'; break;
				case 'T': aa = 'I'; break;
			}
			break;
		}
		break;
	case 'C':
		switch (codon[1])
		{
		case 'A':
			switch (codon[2])
			{
				case 'A': aa = 'Q'; break;
				case 'C': aa = 'H'; break;
				case 'G': aa = 'Q'; break;
				case 'T': aa = 'H'; break;
			}
			break;
		case 'C':
			switch (codon[2])
			{
				case 'A': aa = 'P'; break;
				case 'C': aa = 'P'; break;
				case 'G': aa = 'P'; break;
				case 'T': aa = 'P'; break;
			}
			break;
		case 'G':
			switch (codon[2])
			{
				case 'A': aa = 'R'; break;
				case 'C': aa = 'R'; break;
				case 'G': aa = 'R'; break;
				case 'T': aa = 'R'; break;
			}
			break;
		case 'T':
			switch (codon[2])
			{
				case 'A': aa = 'L'; break;
				case 'C': aa = 'L'; break;
				case 'G': aa = 'L'; break;
				case 'T': aa = 'L'; break;
			}
			break;
		}
		break;
	case 'G':
		switch (codon[1])
		{
		case 'A':
			switch (codon[2])
			{
				case 'A': aa = 'E'; break;
				case 'C': aa = 'D'; break;
				case 'G': aa = 'E'; break;
				case 'T': aa = 'D'; break;
			}
			break;
		case 'C':
			switch (codon[2])
			{
				case 'A': aa = 'A'; break;
				case 'C': aa = 'A'; break;
				case 'G': aa = 'A'; break;
				case 'T': aa = 'A'; break;
			}
			break;
		case 'G':
			switch (codon[2])
			{
				case 'A': aa = 'G'; break;
				case 'C': aa = 'G'; break;
				case 'G': aa = 'G'; break;
				case 'T': aa = 'G'; break;
			}
			break;
		case 'T':
			switch (codon[2])
			{
				case 'A': aa = 'V'; break;
				case 'C': aa = 'V'; break;
				case 'G': aa = 'V'; break;
				case 'T': aa = 'V'; break;
			}
			break;
		}
		break;
	case 'T':
		switch (codon[1])
		{
		case 'A':
			switch (codon[2])
			{
				case 'A': aa = '*'; break;
				case 'C': aa = 'Y'; break;
				case 'G': aa = '*'; break;
				case 'T': aa = 'Y'; break;
			}
			break;
		case 'C':
			switch (codon[2])
			{
				case 'A': aa = 'S'; break;
				case 'C': aa = 'S'; break;
				case 'G': aa = 'S'; break;
				case 'T': aa = 'S'; break;
			}
			break;
		case 'G':
			switch (codon[2])
			{
				case 'A': aa = '*'; break;
				case 'C': aa = 'C'; break;
				case 'G': aa = 'W'; break;
				case 'T': aa = 'C'; break;
			}
			break;
		case 'T':
			switch (codon[2])
			{
				case 'A': aa = 'L'; break;
				case 'C': aa = 'F'; break;
				case 'G': aa = 'L'; break;
				case 'T': aa = 'F'; break;
			}
			break;
		}
		break;
	}
	
	return (aa == '*') ? 0 : aa;
}
