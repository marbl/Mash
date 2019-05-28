// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandTaxScreen.h"
#include "CommandDistance.h" // for pvalue
#include "Sketch.h"
#include "kseq.h"
#include "taxdb.hpp"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"
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

using std::ifstream;
using std::stringstream;

namespace mash {

inline bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}


CommandTaxScreen::CommandTaxScreen()
: Command()
{
	name = "taxscreen";
	summary = "Create Kraken-style taxonomic report based on mash screen.";
	description = "Create Kraken-style taxonomic report based on how well query sequences are contained within a pool of sequences. The queries must be formatted as a single Mash sketch file (.msh), created with the `mash sketch` command. The <pool> files can be contigs or reads, in fasta or fastq, gzipped or not, and \"-\" can be given for <pool> to read from standard input. The <pool> sequences are assumed to be nucleotides, and will be 6-frame translated if the <queries> are amino acids. The output fields are [total percent of hashes, number of contained hashes in the clade, number of contained hashes in the taxon, total number of hashes in the clade, total number of hashes in the taxon, rank, taxonomy ID, padded name].";
    argumentString = "<queries>.msh <pool> [<pool>] ...";

	useOption("help");
	useOption("threads");
//	useOption("minCov");
//    addOption("saturation", Option(Option::Boolean, "s", "", "Include saturation curve in output. Each line will have an additional field representing the absolute number of k-mers seen at each Jaccard increase, formatted as a comma-separated list.", ""));
    addOption("identity", Option(Option::Number, "i", "Output", "Minimum identity to report. Inclusive unless set to zero, in which case only identities greater than zero (i.e. with at least one shared hash) will be reported. Set to -1 to output everything.", "0", -1., 1.));
    addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
	addOption("mapping-file", Option(Option::String, "m", "", "Mapping file from reference name to taxonomy ID", ""));
	addOption("taxonomy-dir", Option(Option::String, "t", "", "Directory containing NCBI taxonomy dump", "."));
}

int CommandTaxScreen::run() const
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

	bool sat = false;//options.at("saturation").active;

    double pValueMax = options.at("pvalue").getArgumentAsNumber();
    double identityMin = options.at("identity").getArgumentAsNumber();
    string taxonomyDir = options.at("taxonomy-dir").argument;
    string mappingFileName = options.at("mapping-file").argument;

    vector<string> refArgVector;
    refArgVector.push_back(arguments[0]);

	Sketch sketch;
    Sketch::Parameters parameters;

    sketch.initFromFiles(refArgVector, parameters);

    string alphabet;
    sketch.getAlphabetAsString(alphabet);
    setAlphabetFromString(parameters, alphabet.c_str());

	parameters.parallelism = options.at("threads").getArgumentAsNumber();
	parameters.kmerSize = sketch.getKmerSize();
	parameters.noncanonical = sketch.getNoncanonical();
	parameters.use64 = sketch.getUse64();
	parameters.preserveCase = sketch.getPreserveCase();
	parameters.seed = sketch.getHashSeed();
	parameters.minHashesPerWindow = sketch.getMinHashesPerWindow();

	HashTable hashTable;
	unordered_map<uint64_t, std::atomic<uint32_t>> hashCounts;
	unordered_map<uint64_t, TaxID> hashTaxIDs;
	unordered_map<uint64_t, list<uint32_t> > saturationByIndex;

	string namesDumpFile = taxonomyDir + "/names.dmp";
	string nodesDumpFile = taxonomyDir + "/nodes.dmp";
	if (!file_exists(namesDumpFile) || !file_exists(nodesDumpFile)) {
		cerr << "Could not find a file names.dmp or nodes.dmp in directory " << taxonomyDir << "\n" 
		     << " To download the required taxonomy files into the current directory, use the following commands:\n"
			 << "   wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\n"
			 << "   tar xvvf taxdump.tar.gz\n"
			 << endl;
		exit(1);

	}
	cerr << "Loading taxonomy files ..." << endl;
	TaxDB taxdb(namesDumpFile, nodesDumpFile); 

	cerr << "Reading mapping file ..." << endl;
	// Read mapping file - not implemented in favor of specifying the taxonomy ID in the reference comment
	vector<TaxID> referenceTaxIDs(sketch.getReferenceCount(), 0);
	if (mappingFileName != "") {
	  	std::ifstream mappingFile(mappingFileName);
	  	if (!mappingFile.is_open())
	    	throw std::runtime_error("unable to open mapping file");

		string referenceID;
		TaxID taxID;
		unordered_map<string, TaxID> refTaxMap;
	  	while (mappingFile >> taxID) 
		{
			mappingFile.ignore(1);
			getline(mappingFile, referenceID, '\n');
			refTaxMap.emplace(referenceID, taxID);
	  	}
		for ( int i = 0; i < sketch.getReferenceCount(); i ++ )
		{
			auto const it = refTaxMap.find(sketch.getReference(i).name);
			if (it == refTaxMap.end()) {
				// No warning? Could still be mapped based on comment
				//cerr << "Could not find taxID for reference " << sketch.getReference(i).name << endl;
			} else {
				referenceTaxIDs[i] = it->second;
			}
		}	
	} 
	for ( int i = 0; i < sketch.getReferenceCount(); i ++ )
	{
		string word;
		TaxID taxID = referenceTaxIDs[i];
		if (taxID == 0) 
		{
			stringstream comment_stream(sketch.getReference(i).comment);
			while (comment_stream >> word) {
				if (word == "taxid") {
					comment_stream >> taxID;
				}
			}
		}
		if (taxID == 0) {
			cerr << "Could not find taxID for reference " << sketch.getReference(i).name << " in comment field or mapping file!" << endl;
		} else {
			//cerr << "Got taxID " << taxID << " for reference " << sketch.getReference(i).name << endl;
			referenceTaxIDs[i] = taxID;
		}
	}
	

	cerr << "Loading " << arguments[0] << "..." << endl;

	// for each reference
	for ( int i = 0; i < sketch.getReferenceCount(); i++ )
	{
		const HashList & hashes = sketch.getReference(i).hashesSorted;

        // for each hash a reference has
		for ( int j = 0; j < hashes.size(); j++ )
		{
			uint64_t hash = hashes.get64() ? hashes.at(j).hash64 : hashes.at(j).hash32;

			if ( hashTable.count(hash) == 0 )
			{
			    // records the counts for each hash
				hashCounts[hash] = 0;
			}

			// save the reference ID for the hash
			hashTable[hash].insert(i);
		}
	}

	cerr << "   " << hashTable.size() << " distinct hashes." << endl;

	unordered_set<MinHashHeap *> minHashHeaps;

	bool trans = (alphabet == alphabetProtein);

/*	if ( ! trans )
	{
		if ( alphabet != alphabetNucleotide )
		{
			cerr << "ERROR: <query> sketch must have nucleotide or amino acid alphabet" << endl;
			exit(1);
		}

		if ( sketch.getNoncanonical() )
		{
			cerr << "ERROR: nucleotide <query> sketch must be canonical" << endl;
			exit(1);
		}
	}
*/

	int queryCount = arguments.size() - 1;
	cerr << (trans ? "Translating from " : "Streaming from ");

	if ( queryCount == 1 )
	{
		cerr << arguments[1];
	}
	else
	{
		cerr << queryCount << " inputs";
	}

	cerr << "..." << endl;

	int kmerSize = parameters.kmerSize;
	int minCov = 1;//options.at("minCov").getArgumentAsNumber();

	ThreadPool<CommandScreen::HashInput, CommandScreen::HashOutput> threadPool(hashSequence, parameters.parallelism);

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

			if ( fps[f - 1] == 0 )
			{
				cerr << "ERROR: could not open " << arguments[f] << endl;
				exit(1);
			}
		}

		kseqs.push_back(kseq_init(fps[f - 1]));
	}

	// perform round-robin, closing files as they end

	int l;
	uint64_t count = 0;
	//uint64_t kmersTotal = 0;
	uint64_t chunkSize = 1 << 20;
	string input;
	input.reserve(chunkSize);
	list<kseq_t *>::iterator it = kseqs.begin();
	//
	while ( true )
	{
		if ( kseqs.begin() == kseqs.end() )
		{
			l = 0;
		}
		else
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
				if ( it == kseqs.end() )
				{
					it = kseqs.begin();
				}
				//continue;
			}
		}

		if ( input.length() + (l >= kmerSize ? l + 1 : 0) > chunkSize || kseqs.begin() == kseqs.end() )
		{
			// chunk big enough or at the end; time to flush

			// buffer this out since kseq will overwrite (deleted by HashInput destructor)
			//
			char * seqCopy = new char[input.length()];
			//
			memcpy(seqCopy, input.c_str(), input.length());

			if ( minHashHeaps.begin() == minHashHeaps.end() )
			{
				minHashHeaps.emplace(new MinHashHeap(sketch.getUse64(), sketch.getMinHashesPerWindow()));
			}

			threadPool.runWhenThreadAvailable(new CommandScreen::HashInput(hashCounts, *minHashHeaps.begin(), seqCopy, input.length(), parameters, trans));

			input = "";

			minHashHeaps.erase(minHashHeaps.begin());

			while ( threadPool.outputAvailable() )
			{
				useThreadOutput(threadPool.popOutputWhenAvailable(), minHashHeaps);
			}
		}

		if ( kseqs.begin() == kseqs.end() )
		{
			break;
		}

		count++;

		if ( l >= kmerSize )
		{
			input.append(1, '*');
			input.append((*it)->seq.s, l);
		}

		it++;

		if ( it == kseqs.end() )
		{
			it = kseqs.begin();
		}
	}

	if (  l != -1 )
	{
		cerr << "\nERROR: reading inputs" << endl;
		exit(1);
	}

	while ( threadPool.running() )
	{
		useThreadOutput(threadPool.popOutputWhenAvailable(), minHashHeaps);
	}

	for ( int i = 0; i < queryCount; i++ )
	{
		gzclose(fps[i]);
	}

	MinHashHeap minHashHeap(sketch.getUse64(), sketch.getMinHashesPerWindow());

	for ( unordered_set<MinHashHeap *>::const_iterator i = minHashHeaps.begin(); i != minHashHeaps.end(); i++ )
	{
		HashList hashList(parameters.use64);

		(*i)->toHashList(hashList);

		for ( int i = 0; i < hashList.size(); i++ )
		{
			minHashHeap.tryInsert(hashList.at(i));
		}

		delete *i;
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
	*/

	uint64_t setSize = minHashHeap.estimateSetSize();
	cerr << "   Estimated distinct" << (trans ? " (translated)" : "") << " k-mers in pool: " << setSize << endl;

	if ( setSize == 0 )
	{
		cerr << "WARNING: no valid k-mers in input." << endl;
		//exit(0);
	}

    // for each hash we can calculate the LCA, and add a count to the LCA at the end
	cerr << "Assigning LCA taxIDs to hashes ..." << endl;

	uint64_t * shared = new uint64_t[sketch.getReferenceCount()]; // number of hashes shared with the query
	vector<uint64_t> * depths = new vector<uint64_t>[sketch.getReferenceCount()]; // how many other references share each hash of a reference?
	memset(shared, 0, sizeof(uint64_t) * sketch.getReferenceCount());
	unordered_map<TaxID, TaxCounts> counts;
	unordered_set<TaxID> allTaxIDs;

	for ( unordered_map<uint64_t, std::atomic<uint32_t> >::const_iterator i = hashCounts.begin(); i != hashCounts.end(); i++ )
	{
		// indices of all the references - map them to taxonomy IDs
		const unordered_set<uint64_t> & indeces = hashTable.at(i->first);

		TaxID taxID = 0;
		for ( unordered_set<uint64_t>::const_iterator k = indeces.begin(); k != indeces.end(); k++ )
		{
			taxID = taxdb.getLowestCommonAncestor(referenceTaxIDs[*k], taxID);
			shared[*k]++;
			depths[*k].push_back(i->second);

			if ( sat )
			{
				saturationByIndex[*k].push_back(0);// TODO kmersTotal);
			}
		}
		//hashTaxIDs.insert(i->first, taxID);
		hashTaxIDs[i->first] = taxID;
		counts[taxID].taxHashCount += 1;
		if ( i->second >= minCov )
		{
			counts[taxID].taxCount += 1;
		allTaxIDs.insert(taxID);
		}
	}

	// Sum up the clade counts and populate the children vectors
	uint64_t totalCount = 0;
	uint64_t totalHashCount = 0;
	for ( unordered_map<TaxID, TaxCounts>::iterator it = counts.begin(); it != counts.end(); ++it ) 
	{
		uint64_t hashCount = it->second.taxHashCount;
		totalHashCount += hashCount;

		uint64_t count = it->second.taxCount;
		totalCount += count;

		TaxID taxID = it->first;
		TaxEntry const * taxon = taxdb.getEntry(taxID);
		while (taxon != NULL) {
			counts[taxon->taxID].cladeCount += count;
			counts[taxon->taxID].cladeHashCount += hashCount;
			if (taxon->parent != NULL) {
				vector<TaxID>& children = counts[taxon->parent->taxID].children;
				auto pc_it = lower_bound(children.begin(),
				                         children.end(),
									     taxon->taxID);
				if (pc_it == children.end() || *pc_it != taxon->taxID) {
					children.insert(pc_it, taxon->taxID);
				}
				taxon = taxon->parent;
			} else {
				break;
			}
		}
	}

	cerr << "Writing output..." << endl;

	taxdb.writeReport(stdout, counts, totalCount, totalHashCount);

	delete [] shared;

	return 0;
}





} // namespace mash
