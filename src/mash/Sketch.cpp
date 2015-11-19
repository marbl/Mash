// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "Sketch.h"
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fcntl.h>
#include <map>
#include "kseq.h"
#include "MurmurHash3.h"
#include <assert.h>
#include <queue>
#include <deque>
#include <set>
#include "Command.h" // TEMP for column printing
#include <sys/stat.h>
#include <capnp/message.h>
#include <capnp/serialize.h>

#define SET_BINARY_MODE(file)
#define CHUNK 16384
KSEQ_INIT(gzFile, gzread)

using namespace std;

typedef map < Sketch::hash_t, vector<Sketch::PositionHash> > LociByHash_map;

const vector<Sketch::Locus> & Sketch::getLociByHash(Sketch::hash_t hash) const
{
    return lociByHash.at(hash);
}

int Sketch::getMinKmerSize(uint64_t reference) const
{
	return ceil(log(references[reference].length * (1 - parameters.warning) / parameters.warning) / log(parameters.protein ? 20 : 4));
}

double Sketch::getRandomKmerChance(uint64_t reference) const
{
	return 1. / (pow(parameters.protein ? 20 : 4, parameters.kmerSize) / references[reference].length + 1.);
}

uint64_t Sketch::getReferenceIndex(string id) const
{
    if ( referenceIndecesById.count(id) == 1 )
    {
        return referenceIndecesById.at(id);
    }
    else
    {
        return -1;
    }
}

void Sketch::initFromBase(const std::string & fileSeq, bool windowed)
{
    file = fileSeq + (windowed ? suffixSketchWindowed : suffixSketch);
    initFromCapnp(file.c_str());
}

int Sketch::initFromCapnp(const char * file, bool headerOnly, bool append)
{
    int fd = open(file, O_RDONLY);
    
    if ( fd < 0 )
    {
        cerr << "ERROR: could not open \"" << file << "\" for reading." << endl;
        exit(1);
    }
    
    capnp::ReaderOptions readerOptions;
    
    readerOptions.traversalLimitInWords = 1000000000000;
    readerOptions.nestingLimit = 1000000;
    
    capnp::StreamFdMessageReader * message = new capnp::StreamFdMessageReader(fd, readerOptions);
    capnp::MinHash::Reader reader = message->getRoot<capnp::MinHash>();
    
    parameters.kmerSize = reader.getKmerSize();
    parameters.error = reader.getError();
    parameters.minHashesPerWindow = reader.getMinHashesPerWindow();
    parameters.windowSize = reader.getWindowSize();
    parameters.concatenated = reader.getConcatenated();
    parameters.noncanonical = reader.getNoncanonical();
    
    if ( headerOnly )
    {
        close(fd);
        
        try
        {
            delete message;
        }
        catch (exception e) {}
        
        return 0;
    }
    
    capnp::MinHash::ReferenceList::Reader referenceListReader = reader.getReferenceList();
    
    capnp::List<capnp::MinHash::ReferenceList::Reference>::Reader referencesReader = referenceListReader.getReferences();
    
    uint64_t referencesOffset = append ? references.size() : 0;
    
    references.resize(referencesOffset + referencesReader.size());
    
    for ( uint64_t i = 0; i < referencesReader.size(); i++ )
    {
        capnp::MinHash::ReferenceList::Reference::Reader referenceReader = referencesReader[i];
        
        Sketch::Reference & reference = references[referencesOffset + i];
        
        reference.name = referenceReader.getName();
        reference.comment = referenceReader.getComment();
        
        if ( referenceReader.getLength64() )
        {
        	reference.length = referenceReader.getLength64();
        }
        else
        {
	        reference.length = referenceReader.getLength();
	    }
        
        reference.hashesSorted.setUse64(parameters.kmerSize > 16);
        
        if ( parameters.kmerSize > 16 )
        {
            capnp::List<uint64_t>::Reader hashesReader = referenceReader.getHashes64();
        
            reference.hashesSorted.resize(hashesReader.size());
        
            for ( uint64_t j = 0; j < hashesReader.size(); j++ )
            {
                reference.hashesSorted.set64(j, hashesReader[j]);
            }
        }
        else
        {
            capnp::List<uint32_t>::Reader hashesReader = referenceReader.getHashes32();
        
            reference.hashesSorted.resize(hashesReader.size());
        
            for ( uint64_t j = 0; j < hashesReader.size(); j++ )
            {
                reference.hashesSorted.set32(j, hashesReader[j]);
            }
        }
    }
    
    capnp::MinHash::LocusList::Reader locusListReader = reader.getLocusList();
    capnp::List<capnp::MinHash::LocusList::Locus>::Reader lociReader = locusListReader.getLoci();
    
    positionHashesByReference.resize(references.size());
    
    for ( uint64_t i = 0; i < lociReader.size(); i++ )
    {
        capnp::MinHash::LocusList::Locus::Reader locusReader = lociReader[i];
        //cout << locusReader.getHash() << '\t' << locusReader.getSequence() << '\t' << locusReader.getPosition() << endl;
        positionHashesByReference[locusReader.getSequence() + referencesOffset].push_back(PositionHash(locusReader.getPosition(), locusReader.getHash64()));
    }
    
    /*
    cout << endl << "References:" << endl << endl;
    
    vector< vector<string> > columns(3);
    
    columns[0].push_back("ID");
    columns[1].push_back("Length");
    columns[2].push_back("Name/Comment");
    
    for ( uint64_t i = 0; i < references.size(); i++ )
    {
        columns[0].push_back(to_string(i));
        columns[1].push_back(to_string(references[i].length));
        columns[2].push_back(references[i].name + " " + references[i].comment);
    }
    
    printColumns(columns);
    cout << endl;
    */
    
    /*
    printf("\nCombined hash table:\n");
    
    cout << "   kmer:  " << kmerSize << endl;
    cout << "   comp:  " << compressionFactor << endl << endl;
    
    for ( LociByHash_umap::iterator i = lociByHash.begin(); i != lociByHash.end(); i++ )
    {
        printf("Hash %u:\n", i->first);
        
        for ( int j = 0; j < i->second.size(); j++ )
        {
            printf("   Seq: %d\tPos: %d\n", i->second.at(j).sequence, i->second.at(j).position);
        }
    }
    
    cout << endl;
    */
    close(fd);
    delete message;
    
    createIndex();
    
    return 0;
}

int Sketch::initFromSequence(const vector<string> & files, const Parameters & parametersNew, int verbosity)
{
    parameters = parametersNew;
    
    int l;
    int count = 0;
    
    bool use64 = parameters.kmerSize > 16;
    
    bloom_filter * bloomFilter = 0;
    
    uint64_t kmersTotal;
    uint64_t kmersUsed;
    
    for ( int i = 0; i < files.size(); i++ )
    {
        HashSet minHashes(parameters.kmerSize);
        HashPriorityQueue minHashesQueue(parameters.kmerSize); // only used for non-windowed
        
        FILE * inStream = 0;
        
        if ( files[i] == "-" )
        {
            if ( verbosity > 0 )
            {
                cerr << "Sketching from stdin..." << endl;
            }
            
            inStream = stdin;
        }
        else
        {
            if ( verbosity > 0 )
            {
                cerr << "Sketching " << files[i] << "..." << endl;
            }
            
            inStream = fopen(files[i].c_str(), "r");
            
            if ( inStream == NULL )
            {
                cerr << "ERROR: could not open " << files[i] << " for reading." << endl;
                exit(1);
            }
        }
        
        gzFile fp = gzdopen(fileno(inStream), "r");
        kseq_t *seq = kseq_init(fp);
        
        if ( parameters.concatenated )
        {
            references.resize(references.size() + 1);
            
            if ( files[i] != "-" )
            {
                references[references.size() - 1].name = files[i];
            }
            
            references[references.size() - 1].length = 0;
            references[references.size() - 1].hashesSorted.setUse64(use64);
            
            if ( parameters.bloomFilter )
            {
            	references[references.size() - 1].length = parameters.genomeSize;
                bloom_parameters bloomParams;
                
                bloomParams.projected_element_count = (uint64_t)parameters.genomeSize * 10l; // TODO: error rate based on platform and coverage
                bloomParams.false_positive_probability = parameters.bloomError;
                bloomParams.maximum_size = (uint64_t)parameters.memoryMax * 8l;
                bloomParams.compute_optimal_parameters();
                
                kmersTotal = 0;
                kmersUsed = 0;
                
                if ( i == 0 && verbosity > 0 )
                {
                    //cerr << "   Bloom table size (bytes): " << bloomParams.optimal_parameters.table_size / 8 << endl;
                }
                
                bloomFilter = new bloom_filter(bloomParams);
            }
        }
        
        bool skipped = false;
        
        while ((l = kseq_read(seq)) >= 0)
        {
            if ( l < parameters.kmerSize )
            {
            	skipped = true;
                continue;
            }
            
            if ( parameters.windowed )
            {
                positionHashesByReference.resize(count + 1);
            }
            
            if ( ! parameters.concatenated )
            {
                references.resize(references.size() + 1);
                references[references.size() - 1].hashesSorted.setUse64(use64);
                minHashes.clear();
                minHashesQueue.clear();
            }
            
            if ( verbosity > 0 && parameters.windowed ) cout << '>' << seq->name.s << " (" << l << "nt)" << endl << endl;
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            
            Reference & reference = references[references.size() - 1];
            
            if ( ! parameters.concatenated )
            {
                reference.name = seq->name.s;
                
                if ( seq->comment.l > 0 )
                {
                    reference.comment = seq->comment.s;
                }
                
                reference.length = l;
            }
            else if ( ! parameters.bloomFilter )
            {
                references[references.size() - 1].length += l;
                
                if ( files[i] == "-" && references[references.size() - 1].name == "" )
                {
                    references[references.size() - 1].name = seq->name.s;
                }
            }
            
            if ( parameters.windowed )
            {
                getMinHashPositions(positionHashesByReference[count], seq->seq.s, l, parameters, verbosity);
            }
            else
            {
                addMinHashes(minHashes, minHashesQueue, bloomFilter, seq->seq.s, l, parameters, kmersTotal, kmersUsed);
            }
            
            if ( ! parameters.concatenated )
            {
                if ( ! parameters.windowed )
                {
                    setMinHashesForReference(references.size() - 1, minHashes);
                }
                
                count++;
            }
        }
        
        if (  l != -1 )
        {
        	cerr << "\nERROR: reading " << files[i] << "." << endl;
        	exit(1);
        }
        
        if ( references[0].length == 0 )
        {
        	if ( skipped )
        	{
        		cerr << "\nWARNING: All fasta records in " << files[i] << "were shorter than the k-mer size (" << parameters.kmerSize << ")." << endl;
        	}
        	else
        	{
        		cerr << "\nERROR: Did not find fasta records in \"" << files[i] << "\"." << endl;
        	}
        	
            exit(1);
        }
        
        if ( parameters.concatenated )
        {
            if ( ! parameters.windowed )
            {
                setMinHashesForReference(references.size() - 1, minHashes);
            }
            
            if ( bloomFilter != 0 )
            {
                if ( verbosity > 0 )
                {
                    //cerr << "   " << kmersTotal - kmersUsed << " of " << kmersTotal << " kmers filtered from " << (files[i] == "-" ? "stdin" : files[i]) << endl;
                }
                
                delete bloomFilter;
            }
            
            count++;
        }
        
        kseq_destroy(seq);
        gzclose(fp);
        fclose(inStream);
    }
    /*
    printf("\nCombined hash table:\n\n");
    
    for ( LociByHash_umap::iterator i = lociByHash.begin(); i != lociByHash.end(); i++ )
    {
        printf("Hash %u:\n", i->first);
        
        for ( int j = 0; j < i->second.size(); j++ )
        {
            printf("   Seq: %d\tPos: %d\n", i->second.at(j).sequence, i->second.at(j).position);
        }
    }
    */
    
    createIndex();
    
    return 0;
}

bool Sketch::initHeaderFromBaseIfValid(const std::string & fileSeq, bool windowed)
{
    file = fileSeq + (windowed ? suffixSketchWindowed : suffixSketch);
    
    struct stat fileInfoSeq;
    struct stat fileInfoSketch;
    
    if ( stat(fileSeq.c_str(), &fileInfoSeq) == -1 )
    {
        return false;
    }
    
    if ( stat(file.c_str(), &fileInfoSketch) == -1 )
    {
        return false;
    }
    
    if ( fileInfoSeq.st_mtime > fileInfoSketch.st_mtime )
    {
        return false;
    }
    
    initFromCapnp(file.c_str(), true);
    return true;
}

void Sketch::warnKmerSize(uint64_t lengthMax, const std::string & lengthMaxName, double randomChance, int kMin, int warningCount) const
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
		<< " is required.\n\n";
}

bool Sketch::writeToFile() const
{
    return writeToCapnp(file.c_str()) == 0;
}

int Sketch::writeToCapnp(const char * file) const
{
    int fd = open(file, O_CREAT | O_WRONLY | O_TRUNC, 0644);
    
    if ( fd < 0 )
    {
        cerr << "ERROR: could not open " << file << " for writing.\n";
        exit(1);
    }
    
    capnp::MallocMessageBuilder message;
    capnp::MinHash::Builder builder = message.initRoot<capnp::MinHash>();
    
    capnp::MinHash::ReferenceList::Builder referenceListBuilder = builder.initReferenceList();
    
    capnp::List<capnp::MinHash::ReferenceList::Reference>::Builder referencesBuilder = referenceListBuilder.initReferences(references.size());
    
    for ( uint64_t i = 0; i < references.size(); i++ )
    {
        capnp::MinHash::ReferenceList::Reference::Builder referenceBuilder = referencesBuilder[i];
        
        referenceBuilder.setName(references[i].name);
        referenceBuilder.setComment(references[i].comment);
        referenceBuilder.setLength64(references[i].length);
        
        if ( references[i].hashesSorted.size() != 0 )
        {
            const HashList & hashes = references[i].hashesSorted;
            
            if ( parameters.kmerSize > 16 )
            {
                capnp::List<uint64_t>::Builder hashes64Builder = referenceBuilder.initHashes64(hashes.size());
            
                uint64_t index = 0;
            
                for ( uint64_t j = 0; j != hashes.size(); j++ )
                {
                    hashes64Builder.set(index, hashes.at(j).hash64);
                    index++;
                }
            }
            else
            {
                capnp::List<uint32_t>::Builder hashes32Builder = referenceBuilder.initHashes32(hashes.size());
            
                uint64_t index = 0;
            
                for ( uint64_t j = 0; j != hashes.size(); j++ )
                {
                    hashes32Builder.set(index, hashes.at(j).hash32);
                    index++;
                }
            }
        }
    }
    
    int locusCount = 0;
    
    for ( int i = 0; i < positionHashesByReference.size(); i++ )
    {
        locusCount += positionHashesByReference.at(i).size();
    }
    
    capnp::MinHash::LocusList::Builder locusListBuilder = builder.initLocusList();
    capnp::List<capnp::MinHash::LocusList::Locus>::Builder lociBuilder = locusListBuilder.initLoci(locusCount);
    
    int locusIndex = 0;
    
    for ( int i = 0; i < positionHashesByReference.size(); i++ )
    {
        for ( int j = 0; j < positionHashesByReference.at(i).size(); j++ )
        {
            capnp::MinHash::LocusList::Locus::Builder locusBuilder = lociBuilder[locusIndex];
            locusIndex++;
            
            locusBuilder.setSequence(i);
            locusBuilder.setPosition(positionHashesByReference.at(i).at(j).position);
            locusBuilder.setHash64(positionHashesByReference.at(i).at(j).hash);
        }
    }
    
    builder.setKmerSize(parameters.kmerSize);
    builder.setError(parameters.error);
    builder.setMinHashesPerWindow(parameters.minHashesPerWindow);
    builder.setWindowSize(parameters.windowSize);
    builder.setConcatenated(parameters.concatenated);
    builder.setNoncanonical(parameters.noncanonical);
    
    writeMessageToFd(fd, message);
    close(fd);
    
    return 0;
}

void Sketch::createIndex()
{
    for ( int i = 0; i < references.size(); i++ )
    {
        referenceIndecesById[references[i].name] = i;
    }
    
    for ( int i = 0; i < positionHashesByReference.size(); i++ )
    {
        for ( int j = 0; j < positionHashesByReference.at(i).size(); j++ )
        {
            const PositionHash & positionHash = positionHashesByReference.at(i).at(j);
            
            lociByHash[positionHash.hash].push_back(Locus(i, positionHash.position));
        }
    }
    
    kmerSpace = pow(4, parameters.kmerSize); // TODO: alphabet?
}

void Sketch::setMinHashesForReference(uint64_t referenceIndex, const HashSet & hashes)
{
    HashList & hashList = references[referenceIndex].hashesSorted;
    hashList.clear();
    hashes.toHashList(hashList);
    hashList.sort();
}

void addMinHashes(HashSet & minHashes, HashPriorityQueue & minHashesQueue, bloom_filter * bloomFilter, char * seq, uint64_t length, const Sketch::Parameters & parameters, uint64_t & kmersTotal, uint64_t & kmersUsed)
{
    int kmerSize = parameters.kmerSize;
    uint64_t mins = parameters.minHashesPerWindow;
    bool noncanonical = parameters.noncanonical;
    
    // Determine the 'mins' smallest hashes, including those already provided
    // (potentially replacing them). This allows min-hash sets across multiple
    // sequences to be determined.
    
    // uppercase TODO: alphabets?
    //
    for ( uint64_t i = 0; i < length; i++ )
    {
        if ( seq[i] > 90 )
        {
            seq[i] -= 32;
        }
    }
    
    bool use64 = kmerSize > 16;
    
    char * seqRev = new char[length];
    
    if ( ! noncanonical )
    {
        reverseComplement(seq, seqRev, length);
    }
    
    for ( uint64_t i = 0; i < length - kmerSize + 1; i++ )
    {
        bool useRevComp = false;
        bool debug = false;
        
        if ( ! noncanonical )
        {
            // repeatedly skip kmers with bad characters
            
            bool bad = false;
            
            for ( uint64_t j = i; j < i + kmerSize && i + kmerSize <= length; j++ )
            {
                char c = seq[j];
            
                if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' )
                {
                    i = j; // skip to past the bad character
                    bad = true;
                    break;
                }
            }
            
            if ( bad )
            {
                continue;
            }
        
            if ( i + kmerSize > length )
            {
                // skipped to end
                break;
            }
            
            useRevComp = true;
            bool prefixEqual = true;
        
            if ( debug ) {for ( uint64_t j = i; j < i + kmerSize; j++ ) { cout << *(seq + j); } cout << endl;}
        
            for ( uint64_t j = 0; j < kmerSize; j++ )
            {
                char base = seq[i + j];
                char baseMinus = seqRev[length - i - kmerSize + j];
            
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
        
        const char * kmer = useRevComp ? seqRev + length - i - kmerSize : seq + i;
        bool filter = false;
        
        hash_u hash = getHash(useRevComp ? seqRev + length - i - kmerSize : seq + i, kmerSize);
        
        if ( debug ) cout << endl;
        
        if
        (
            (
                minHashesQueue.size() < mins ||
                hashLessThan(hash, minHashesQueue.top(), use64)
            )
            && ! minHashes.contains(hash)
        )
        {
            if ( bloomFilter )
            {
                kmersTotal++;
                std::string kmerString(kmer, kmerSize);
            
                //cout << kmerString;
            
                if ( ! bloomFilter->contains(kmerString) )
                {
                    filter = true;
                    //cout << " (filtered)";
                }
            
                //cout << endl;
            
                bloomFilter->insert(kmerString);
            
                if ( filter )
                {
                    continue;
                }
            
                kmersUsed++;
            }
        
            minHashes.insert(hash);
            minHashesQueue.push(hash);
            
            if ( minHashesQueue.size() > mins )
            {
                minHashes.erase(minHashesQueue.top());
                minHashesQueue.pop();
            }
        }
    }
    
    if ( ! noncanonical )
    {
        delete [] seqRev;
    }
}

void getMinHashPositions(vector<Sketch::PositionHash> & positionHashes, char * seq, uint32_t length, const Sketch::Parameters & parameters, int verbosity)
{
    // Find positions whose hashes are min-hashes in any window of a sequence
    
    int kmerSize = parameters.kmerSize;
    int mins = parameters.minHashesPerWindow;
    int windowSize = parameters.windowSize;
    
    int nextValidKmer = 0;
    
    if ( windowSize > length - kmerSize + 1 )
    {
        windowSize = length - kmerSize + 1;
    }
    
    if ( verbosity > 1 ) cout << seq << endl << endl;
    
    // Associate positions with flags so they can be marked as min-hashes
    // at any point while the window is moved across them
    //
    struct CandidateLocus
    {
        CandidateLocus(int positionNew)
            :
            position(positionNew),
            isMinmer(false)
            {}
        
        int position;
        bool isMinmer;
    };
    
    // All potential min-hash loci in the current window organized by their
    // hashes so repeats can be grouped and so the sorted keys can be used to
    // keep track of the current h bottom hashes. A deque is used here (rather
    // than a standard queue) for each list of candidate loci for the sake of
    // debug output; the performance difference seems to be negligible.
    //
    map<Sketch::hash_t, deque<CandidateLocus>> candidatesByHash;
    
    // Keep references to the candidate loci in the map in the order of the
    // current window, allowing them to be popped off in the correct order as
    // the window is incremented.
    //
    queue<map<Sketch::hash_t, deque<CandidateLocus>>::iterator> windowQueue;
    
    // Keep a reference to the "hth" min-hash to determine efficiently whether
    // new hashes are min-hashes. It must be decremented when a hash is inserted
    // before it.
    //
    map<Sketch::hash_t, deque<CandidateLocus>>::iterator maxMinmer = candidatesByHash.end();
    
    // A reference to the list of potential candidates that has just been pushed
    // onto the back of the rolling window. During the loop, it will be assigned
    // to either an existing list (if the kmer is repeated), a new list, or a
    // dummy iterator (for invalid kmers).
    //
    map<Sketch::hash_t, deque<CandidateLocus>>::iterator newCandidates;
    
    int unique = 0;
    
    for ( int i = 0; i < length - kmerSize + 1; i++ )
    {
        // Increment the next valid kmer if needed. Invalid kmers must still be
        // processed to keep the queue filled, but will be associated with a
        // dummy iterator. (Currently disabled to allow all kmers; see below)
        //
        if ( i >= nextValidKmer )
        {
            for ( int j = i; j < i + kmerSize; j++ )
            {
                char c = seq[j];
                
                if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' )
                {
                    // Uncomment to skip invalid kmers
                    //
                    //nextValidKmer = j + 1;
                    
                    break;
                }
            }
        }
        
        if ( i < nextValidKmer && verbosity > 1 )
        {
            cout << "  [";
        
            for ( int j = i; j < i + kmerSize; j++ )
            {
                cout << seq[j];
            }
            
            cout << "]" << endl;
        }
        
        if ( i >= nextValidKmer )
        {
            Sketch::hash_t hash = getHash(seq + i, kmerSize).hash64; // TODO: dynamic
            
            if ( verbosity > 1 )
            {
                cout << "   ";
            
                for ( int j = i; j < i + kmerSize; j++ )
                {
                    cout << seq[j]; 
                }
            
                cout << "   " << i << '\t' << hash << endl;
            }
            
            // Get the list of candidate loci for the current hash (if it is a
            // repeat) or insert a new list.
            //
            pair<map<Sketch::hash_t, deque<CandidateLocus>>::iterator, bool> inserted =
                candidatesByHash.insert(pair<Sketch::hash_t, deque<CandidateLocus>>(hash, deque<CandidateLocus>()));
            newCandidates = inserted.first;
            
            // Add the new candidate locus to the list
            //
            newCandidates->second.push_back(CandidateLocus(i));
            
            if
            (
                inserted.second && // inserted; decrement maxMinmer if...
                (
                    (
                        // ...just reached number of mins
                        
                        maxMinmer == candidatesByHash.end() &&
                        candidatesByHash.size() == mins
                    ) ||
                    (
                        // ...inserted before maxMinmer
                        
                        maxMinmer != candidatesByHash.end() &&
                        newCandidates->first < maxMinmer->first
                    )
                )
            )
            {
                maxMinmer--;
                
                if ( i >= windowSize )
                {
                    unique++;
                }
            }
        }
        else
        {
            // Invalid kmer; use a dummy iterator to pad the queue
            
            newCandidates = candidatesByHash.end();
        }
        
        // Push the new reference to the list of candidate loci for the new kmer
        // on the back of the window to roll it.
        //
        windowQueue.push(newCandidates);
        
        // A reference to the front of the window, to be popped off if the
        // window has grown to full size. This reference can be a dummy if the
        // window is not full size or if the front of the window is a dummy
        // representing an invalid kmer.
        //
        map<Sketch::hash_t, deque<CandidateLocus>>::iterator windowFront = candidatesByHash.end();
        
        if ( windowQueue.size() > windowSize )
        {
            windowFront = windowQueue.front();
            windowQueue.pop();
            
            if ( verbosity > 1 ) cout << "   \tPOP: " << windowFront->first << endl;
        }
        
        if ( windowFront != candidatesByHash.end() )
        {
            deque<CandidateLocus> & frontCandidates = windowFront->second;
            
            if ( frontCandidates.front().isMinmer )
            {
                if ( verbosity > 1 ) cout << "   \t   minmer: " << frontCandidates.front().position << '\t' << windowFront->first << endl;
                positionHashes.push_back(Sketch::PositionHash(frontCandidates.front().position, windowFront->first));
            }
            
            if ( frontCandidates.size() > 1 )
            {
                frontCandidates.pop_front();
                
                // Since this is a repeated hash, only the locus in the front of
                // the list was considered min-hash loci. Check if the new front
                // will become a min-hash so it can be flagged.
                //
                if ( maxMinmer == candidatesByHash.end() || ( i >= windowSize && windowFront->first <= maxMinmer->first) )
                {
                    frontCandidates.front().isMinmer = true;
                }
            }
            else
            {
                // The list for this hash is no longer needed; destroy it,
                // repositioning the reference to the hth min-hash if
                // necessary.
                
                if ( maxMinmer != candidatesByHash.end() && windowFront->first <= maxMinmer->first )
                {
                    maxMinmer++;
                    
                    if ( maxMinmer != candidatesByHash.end() )
                    {
                        maxMinmer->second.front().isMinmer = true;
                    }
                    
                    unique++;
                }
            
                candidatesByHash.erase(windowFront);
            }
        }
        
        if ( i == windowSize - 1 )
        {
            // first complete window; mark min-hashes
            
            for ( map<Sketch::hash_t, deque<CandidateLocus>>::iterator j = candidatesByHash.begin(); j != maxMinmer; j++ )
            {
                j->second.front().isMinmer = true;
            }
            
            if ( maxMinmer != candidatesByHash.end() )
            {
                maxMinmer->second.front().isMinmer = true;
            }
            
            unique++;
        }
        
        // Mark the candidate that was pushed on the back of the window queue
        // earlier as a min-hash if necessary
        //
        if ( newCandidates != candidatesByHash.end() && i >= windowSize && (maxMinmer == candidatesByHash.end() || newCandidates->first <= maxMinmer->first) )
        {
            newCandidates->second.front().isMinmer = true;
        }
        
        if ( verbosity > 1 )
        {
            for ( map<Sketch::hash_t, deque<CandidateLocus>>::iterator j = candidatesByHash.begin(); j != candidatesByHash.end(); j++ )
            {
                cout << "   \t" << j->first;
                
                if ( j == maxMinmer )
                {
                     cout << "*";
                }
                
                for ( deque<CandidateLocus>::iterator k = j->second.begin(); k != j->second.end(); k++ )
                {
                    cout << '\t' << k->position;
                    
                    if ( k->isMinmer )
                    {
                        cout << '!';
                    }
                }
                
                cout << endl;
            }
        }
    }
    
    // finalize remaining min-hashes from the last window
    //
    while ( windowQueue.size() > 0 )
    {
        map<Sketch::hash_t, deque<CandidateLocus>>::iterator windowFront = windowQueue.front();
        windowQueue.pop();
        
        if ( windowFront != candidatesByHash.end() )
        {
            deque<CandidateLocus> & frontCandidates = windowFront->second;
            
            if ( frontCandidates.size() > 0 )
            {
                if ( frontCandidates.front().isMinmer )
                {
                    if ( verbosity > 1 ) cout << "   \t   minmer:" << frontCandidates.front().position << '\t' << windowFront->first << endl;
                    positionHashes.push_back(Sketch::PositionHash(frontCandidates.front().position, windowFront->first));
                }
                
                frontCandidates.pop_front();
            }
        }
    }
    
    if ( verbosity > 1 )
    {
        cout << endl << "Minmers:" << endl;
    
        for ( int i = 0; i < positionHashes.size(); i++ )
        {
            cout << "   " << positionHashes.at(i).position << '\t' << positionHashes.at(i).hash << endl;
        }
        
        cout << endl;
    }
    
    if ( verbosity > 0 ) cout << "   " << positionHashes.size() << " minmers across " << length - windowSize - kmerSize + 2 << " windows (" << unique << " windows with distinct minmer sets)." << endl << endl;
}

bool hasSuffix(string const & whole, string const & suffix)
{
    if (whole.length() >= suffix.length())
    {
        return 0 == whole.compare(whole.length() - suffix.length(), suffix.length(), suffix);
    }
    
    return false;
}

void reverseComplement(const char * src, char * dest, int length)
{
    for ( int i = 0; i < length; i++ )
    {
        char base = src[i];
        
        switch ( base )
        {
            case 'A': base = 'T'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
            case 'T': base = 'A'; break;
            default: break;
        }
        
        dest[length - i - 1] = base;
    }
}

// The following functions are adapted from http://www.zlib.net/zpipe.c


/* Compress from file source to file dest until EOF on source.
   def() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_STREAM_ERROR if an invalid compression
   level is supplied, Z_VERSION_ERROR if the version of zlib.h and the
   version of the library linked do not match, or Z_ERRNO if there is
   an error reading or writing the files. */
int def(int fdSource, int fdDest, int level)
{
    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        return ret;

    /* compress until end of file */
    do {
        strm.avail_in = read(fdSource, in, CHUNK);
        if (strm.avail_in == -1) {
            (void)deflateEnd(&strm);
            return Z_ERRNO;
        }
        flush = strm.avail_in == 0 ? Z_FINISH : Z_NO_FLUSH;
        strm.next_in = in;

        /* run deflate() on input until output buffer not full, finish
           compression if all of source has been read in */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = deflate(&strm, flush);    /* no bad return value */
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            have = CHUNK - strm.avail_out;
            if (write(fdDest, out, have) != have) {
                (void)deflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);
        assert(strm.avail_in == 0);     /* all input will be used */

        /* done when last data in file processed */
    } while (flush != Z_FINISH);
    assert(ret == Z_STREAM_END);        /* stream will be complete */

    /* clean up and return */
    (void)deflateEnd(&strm);
    return Z_OK;
}

/* Decompress from file source to file dest until stream ends or EOF.
   inf() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_DATA_ERROR if the deflate data is
   invalid or incomplete, Z_VERSION_ERROR if the version of zlib.h and
   the version of the library linked do not match, or Z_ERRNO if there
   is an error reading or writing the files. */
int inf(int fdSource, int fdDest)
{
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return ret;

    /* decompress until deflate stream ends or end of file */
    do {
        strm.avail_in = read(fdSource, in, CHUNK);
        if (strm.avail_in == -1) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            if (write(fdDest, out, have) != have) {
                (void)inflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

/* report a zlib or i/o error */
void zerr(int ret)
{
    fputs("zpipe: ", stderr);
    switch (ret) {
    case Z_ERRNO:
        if (ferror(stdin))
            fputs("error reading stdin\n", stderr);
        if (ferror(stdout))
            fputs("error writing stdout\n", stderr);
        break;
    case Z_STREAM_ERROR:
        fputs("invalid compression level\n", stderr);
        break;
    case Z_DATA_ERROR:
        fputs("invalid or incomplete deflate data\n", stderr);
        break;
    case Z_MEM_ERROR:
        fputs("out of memory\n", stderr);
        break;
    case Z_VERSION_ERROR:
        fputs("zlib version mismatch!\n", stderr);
    }
}
