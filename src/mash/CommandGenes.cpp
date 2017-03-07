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
	name = "within";
	summary = "Determine whether query sequences are within a larger pool of sequences.";
	description = "Determine whether query sequences are within a larger pool of sequences. The targets must be formatted as a Mash sketch file (.msh), created with the `mash sketch` command. The <pool> files can be contigs or reads, in fasta or fastq, gzipped or not. The output fields are [seq-ID-1, seq-ID-2, distance, p-value, shared-hashes].";
    argumentString = "<queries> <pool> [<pool>] ...";
	
	useOption("help");
	useOption("minCov");
	//useSketchOptions();
}

int CommandGenes::run() const
{
	if ( arguments.size() < 2 || options.at("help").active )
	{
		print();
		return 0;
	}
	
	if ( ! hasSuffix(arguments[0], suffixSketch) )
	{
		cerr << "ERROR: " << arguments[0] << " does not look like a sketch (.msh)" << endl;
		exit(1);
	}
	
    vector<string> refArgVector;
    refArgVector.push_back(arguments[0]);
	
	Sketch sketch;
    Sketch::Parameters parameters;
	
    sketch.initFromFiles(refArgVector, parameters);
    
    string alphabet;
    sketch.getAlphabetAsString(alphabet);
    setAlphabetFromString(parameters, alphabet.c_str());
	
	HashTable hashTable;
	unordered_map<uint64_t, uint32_t> hashCounts;
	
	cerr << "Filling table from " << arguments[0] << endl;
	
	for ( int i = 0; i < sketch.getReferenceCount(); i++ )
	{
		const HashList & hashes = sketch.getReference(i).hashesSorted;
		
		for ( int j = 0; j < hashes.size(); j++ )
		{
			uint64_t hash = hashes.get64() ? hashes.at(j).hash64 : hashes.at(j).hash32;
			hashTable[hash].insert(i);
		}
	}
	
	uint64_t * shared = new uint64_t[sketch.getReferenceCount()];
	
	memset(shared, 0, sizeof(uint64_t) * sketch.getReferenceCount());
	
	MinHashHeap minHashHeap(parameters.use64, parameters.minHashesPerWindow, parameters.minCov, parameters.memoryBound);
	
	int queryCount = arguments.size() - 1;
	cerr << "Streaming from " << queryCount << " inputs..." << endl;
	
	bool trans = false;
	bool use64 = sketch.getUse64();
	int kmerSize = sketch.getKmerSize();
	int minCov = options.at("minCov").getArgumentAsNumber();
	bool noncanonical = sketch.getNoncanonical();
	
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
	int l;
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
		
		if ( l < kmerSize ) // too short
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
		
		char * seqRev;
		
		if ( ! noncanonical || trans )
		{
			seqRev = new char[l];
			reverseComplement(seq, seqRev, l);
		}
		
		for ( int i = 0; i < (trans ? 6 : 1); i++ )
		{
			bool useRevComp = false;
			int frame = i % 3;
			bool rev = i > 2;
			
			int lenTrans = (l - frame) / 3;
			
			char * seqTrans;
			
			if ( trans )
			{
				seqTrans = new char[lenTrans];
				translate((rev ? seqRev : seq) + frame, seqTrans, lenTrans);
			}
			
			int64_t lastGood = -1;
			int length = trans ? lenTrans : l;
			
			for ( int j = 0; j < length - kmerSize + 1; j++ )
			{
				while ( lastGood < j + kmerSize - 1 && lastGood < length )
				{
					lastGood++;
					
					if ( trans ? (seqTrans[lastGood] == 0) : (!parameters.alphabet[seq[lastGood]]) )
					{
						j = lastGood + 1;
					}
				}
				
				if ( j > length - kmerSize )
				{
					break;
				}
				
				if ( ! noncanonical )
				{
					bool debug = false;
					useRevComp = true;
					bool prefixEqual = true;
		
					if ( debug ) {for ( uint64_t k = j; k < j + kmerSize; k++ ) { cout << *(seq + k); } cout << endl;}
					
					for ( uint64_t k = 0; k < kmerSize; k++ )
					{
						char base = seq[j + k];
						char baseMinus = seqRev[l - j - kmerSize + k];
			
						if ( debug ) cout << baseMinus;
			
						if ( prefixEqual && baseMinus > base )
						{
							useRevComp = false;
							break;
						}
			
						if ( prefixEqual && baseMinus < base )
						{
							prefixEqual = false;
						}
					}
		
					if ( debug ) cout << endl;
				}
		
				const char * kmer;
				
				if ( trans )
				{
					kmer = seqTrans + i;
				}
				else
				{
					kmer = useRevComp ? seqRev + l - j - kmerSize : seq + j;
				}
				
				hash_u hash = getHash(kmer, kmerSize, use64);
				
				uint64_t key = use64 ? hash.hash64 : hash.hash32;
				
				if ( hashTable.count(key) == 1 )
				{
					hashCounts[key]++;
					
					if ( hashCounts.at(key) == minCov )
					{
						const unordered_set<uint64_t> & indeces = hashTable.at(key);
				
						for ( unordered_set<uint64_t>::const_iterator k = indeces.begin(); k != indeces.end(); k++ )
						{
							shared[*k]++;
						}
					}
				}
			}
			
			if ( trans )
			{
				delete [] seqTrans;
			}
		}
		
		if ( ! sketch.getNoncanonical() || trans )
		{
			delete [] seqRev;
		}
		/*
		addMinHashes(minHashHeap, seq, l, parameters);
		
		if ( parameters.targetCov > 0 && minHashHeap.estimateMultiplicity() >= parameters.targetCov )
		{
			l = -1; // success code
			break;
		}
		*/
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
	/*
	if ( parameters.targetCov != 0 )
	{
		cerr << "Reads required for " << parameters.targetCov << "x coverage: " << count << endl;
	}
	else
	{
		cerr << "Estimated coverage: " << minHashHeap.estimateMultiplicity() << "x" << endl;
	}
	cerr << "Estimated genome size: " << minHashHeap.estimateSetSize() << endl;
	*/
	for ( int i = 0; i < sketch.getReferenceCount(); i++ )
	{
		if ( shared[i] != 0 )
		{
			double identity = estimateIdentity(shared[i], sketch.getReference(i).hashesSorted.size(), kmerSize, sketch.getKmerSpace());
			
			cout << identity << '\t' << shared[i] << '/' << sketch.getReference(i).hashesSorted.size() << '\t' << sketch.getReference(i).name << '\t' << sketch.getReference(i).comment << endl;
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
		distance = -log(jaccard) / kmerSize;
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
