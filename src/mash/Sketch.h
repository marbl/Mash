// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef Sketch_h
#define Sketch_h

#include "mash/capnp/MinHash.capnp.h"
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include "MinHashHeap.h"
#include "ThreadPool.h"

static const char * capnpHeader = "Cap'n Proto";
static const int capnpHeaderLength = strlen(capnpHeader);

static const char * suffixSketch = ".msh";
static const char * suffixSketchWindowed = ".msw";

static const char * alphabetNucleotide = "ACGT";
static const char * alphabetProtein = "ACDEFGHIKLMNPQRSTVWY";

class Sketch
{
public:
    
    typedef uint64_t hash_t;
    
    struct Parameters
    {
        Parameters()
            :
            parallelism(0),
            kmerSize(0),
            alphabetSize(0),
            preserveCase(false),
            use64(false),
            seed(0),
            error(0),
            warning(0),
            minHashesPerWindow(0),
            windowSize(0),
            windowed(false),
            concatenated(false),
            noncanonical(false),
            reads(false),
            memoryBound(0),
            minCov(0),
            targetCov(0),
            genomeSize(0)
        {
        	memset(alphabet, 0, 256);
        }
        
        Parameters(const Parameters & other)
            :
            parallelism(other.parallelism),
            kmerSize(other.kmerSize),
            alphabetSize(other.alphabetSize),
            preserveCase(other.preserveCase),
            use64(other.use64),
            seed(other.seed),
            error(other.error),
            warning(other.warning),
            minHashesPerWindow(other.minHashesPerWindow),
            windowSize(other.windowSize),
            windowed(other.windowed),
            concatenated(other.concatenated),
            noncanonical(other.noncanonical),
            reads(other.reads),
            memoryBound(other.memoryBound),
            minCov(other.minCov),
            targetCov(other.targetCov),
            genomeSize(other.genomeSize)
		{
			memcpy(alphabet, other.alphabet, 256);
		}
        
        int parallelism;
        int kmerSize;
        bool alphabet[256];
        uint32_t alphabetSize;
        bool preserveCase;
        bool use64;
        uint32_t seed;
        double error;
        double warning;
        uint64_t minHashesPerWindow;
        uint64_t windowSize;
        bool windowed;
        bool concatenated;
        bool noncanonical;
        bool reads;
        uint64_t memoryBound;
        uint32_t minCov;
        double targetCov;
        uint64_t genomeSize;
    };
    
    struct PositionHash
    {
        PositionHash(uint32_t positionNew, hash_t hashNew) :
            position(positionNew),
            hash(hashNew)
            {}

        uint32_t position;
        hash_t hash;
    };
    
    struct Locus
    {
        Locus(uint32_t sequenceNew, uint32_t positionNew)
            :
            sequence(sequenceNew),
            position(positionNew)
            {}
        
        uint32_t sequence;
        uint32_t position;
    };
    
    typedef std::unordered_set<hash_t> Hash_set;
    
    struct Reference
    {
        // no sequence for now
        
        std::string name;
        std::string comment;
        uint64_t length;
        HashList hashesSorted;
        std::vector<uint32_t> counts;
    };
    
    struct SketchInput
    {
    	SketchInput(std::string fileNameNew, char * seqNew, uint64_t lengthNew, const std::string & nameNew, const std::string & commentNew, const Sketch::Parameters & parametersNew)
    	:
    	fileName(fileNameNew),
    	seq(seqNew),
    	length(lengthNew),
    	name(nameNew),
    	comment(commentNew),
    	parameters(parametersNew)
    	{}
    	
    	~SketchInput()
    	{
    		if ( seq != 0 )
    		{
	    		delete [] seq;
	    	}
    	}
    	
    	std::string fileName;
    	
    	char * seq;
    	
    	uint64_t length;
    	
    	std::string name;
    	std::string comment;
    	
    	Sketch::Parameters parameters;
    };
    
    struct SketchOutput
    {
    	std::vector<Reference> references;
	    std::vector<std::vector<PositionHash>> positionHashesByReference;
    };
    
    void getAlphabetAsString(std::string & alphabet) const;
    uint32_t getAlphabetSize() const {return parameters.alphabetSize;}
    bool getConcatenated() const {return parameters.concatenated;}
    float getError() const {return parameters.error;}
    int getHashCount() const {return lociByHash.size();}
    uint32_t getHashSeed() const {return parameters.seed;}
    const std::vector<Locus> & getLociByHash(hash_t hash) const;
    float getMinHashesPerWindow() const {return parameters.minHashesPerWindow;}
	int getMinKmerSize(uint64_t reference) const;
	bool getPreserveCase() const {return parameters.preserveCase;}
	double getRandomKmerChance(uint64_t reference) const;
    const Reference & getReference(uint64_t index) const {return references.at(index);}
    uint64_t getReferenceCount() const {return references.size();}
    void getReferenceHistogram(uint64_t index, std::map<uint32_t, uint64_t> & histogram) const;
    uint64_t getReferenceIndex(std::string id) const;
    int getKmerSize() const {return parameters.kmerSize;}
    double getKmerSpace() const {return kmerSpace;}
    bool getUse64() const {return parameters.use64;}
    uint64_t getWindowSize() const {return parameters.windowSize;}
    bool getNoncanonical() const {return parameters.noncanonical;}
    bool hasHashCounts() const {return references.size() > 0 && references.at(0).counts.size() > 0;}
    bool hasLociByHash(hash_t hash) const {return lociByHash.count(hash);}
    int initFromFiles(const std::vector<std::string> & files, const Parameters & parametersNew, int verbosity = 0, bool enforceParameters = false, bool contain = false);
    void initParametersFromCapnp(const char * file);
	bool sketchFileBySequence(FILE * file, ThreadPool<Sketch::SketchInput, Sketch::SketchOutput> * threadPool);
	void useThreadOutput(SketchOutput * output);
    void warnKmerSize(uint64_t lengthMax, const std::string & lengthMaxName, double randomChance, int kMin, int warningCount) const;
    bool writeToFile() const;
    int writeToCapnp(const char * file) const;
    
private:
    
    void createIndex();
    
    std::vector<Reference> references;
    std::unordered_map<std::string, int> referenceIndecesById;
    std::vector<std::vector<PositionHash>> positionHashesByReference;
    std::unordered_map<hash_t, std::vector<Locus>> lociByHash;
    
    Parameters parameters;
    double kmerSpace;
    std::string file;
};

void addMinHashes(MinHashHeap & minHashHeap, char * seq, uint64_t length, const Sketch::Parameters & parameters);
void getMinHashPositions(std::vector<Sketch::PositionHash> & loci, char * seq, uint32_t length, const Sketch::Parameters & parameters, int verbosity = 0);
bool hasSuffix(std::string const & whole, std::string const & suffix);
Sketch::SketchOutput * loadCapnp(Sketch::SketchInput * input);
void reverseComplement(const char * src, char * dest, int length);
void setAlphabetFromString(Sketch::Parameters & parameters, const char * characters);
void setMinHashesForReference(Sketch::Reference & reference, const MinHashHeap & hashes);
Sketch::SketchOutput * sketchFile(Sketch::SketchInput * input);
Sketch::SketchOutput * sketchSequence(Sketch::SketchInput * input);

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);

#endif
