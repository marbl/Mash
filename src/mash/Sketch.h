#ifndef Sketch_h
#define Sketch_h

#include "mash/capnp/MinHash.capnp.h"
#include <unordered_map>
#include <vector>
#include <string>

#include "HashList.h"
#include "HashPriorityQueue.h"
#include "HashSet.h"
#include <stdlib.h> // needed (but not included) by bloom_filter.hpp
#include "bloom_filter.hpp"

static const char * capnpHeader = "Cap'n Proto";
static const int capnpHeaderLength = strlen(capnpHeader);

static const char * suffixSketch = ".msh";
static const char * suffixSketchWindowed = ".msw";

class Sketch
{
public:
    
    typedef uint64_t hash_t;
    
    struct Parameters
    {
        Parameters()
            :
            kmerSize(0),
            error(0),
            minHashesPerWindow(0),
            windowSize(0),
            windowed(false),
            concatenated(false),
            noncanonical(false),
            bloomFilter(false),
            genomeSize(0),
            memoryMax(0),
            bloomError(0)
            {}
        
        Parameters(const Parameters & other)
            :
            kmerSize(other.kmerSize),
            error(other.error),
            minHashesPerWindow(other.minHashesPerWindow),
            windowSize(other.windowSize),
            windowed(other.windowed),
            concatenated(other.concatenated),
            noncanonical(other.noncanonical),
            bloomFilter(other.bloomFilter),
            genomeSize(other.genomeSize),
            memoryMax(other.memoryMax),
            bloomError(other.bloomError)
            {}
        
        int kmerSize;
        float error;
        int minHashesPerWindow;
        int windowSize;
        bool windowed;
        bool concatenated;
        bool noncanonical;
        bool bloomFilter;
        int genomeSize;
        int memoryMax;
        float bloomError;
    };
    
    struct PositionHash
    {
        PositionHash(uint32_t positionNew, uint32_t hashNew) :
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
    
    typedef std::unordered_map < hash_t, std::vector<Sketch::PositionHash> > LociByHash_umap;
    typedef std::unordered_set<hash_t> Hash_set;
    
    struct Reference
    {
        // no sequence for now
        
        std::string name;
        std::string comment;
        uint32_t length;
        HashList hashesSorted;
    };
    
    bool getConcatenated() const {return parameters.concatenated;}
    float getError() const {return parameters.error;}
    int getHashCount() const {return lociByHash.size();}
    const std::vector<Locus> & getLociByHash(hash_t hash) const;
    float getMinHashesPerWindow() const {return parameters.minHashesPerWindow;}
    const Reference & getReference(int index) const {return references.at(index);}
    int getReferenceCount() const {return references.size();}
    int getReferenceIndex(std::string id) const;
    int getKmerSize() const {return parameters.kmerSize;}
    double getKmerSpace() const {return kmerSpace;}
    int getWindowSize() const {return parameters.windowSize;}
    bool getNoncanonical() const {return parameters.noncanonical;}
    bool hasLociByHash(hash_t hash) const {return lociByHash.count(hash);}
    void initFromBase(const std::string & file, bool windowed);
    int initFromCapnp(const char * file, bool headerOnly = false, bool append = false);
    int initFromSequence(const std::vector<std::string> & files, const Parameters & parametersNew, int verbosity = 0);
    bool initHeaderFromBaseIfValid(const std::string & file, bool windowed);
    bool writeToFile() const;
    int writeToCapnp(const char * file) const;
    
private:
    
    void createIndex();
    void setMinHashesForReference(int referenceIndex, const HashSet & hashes);
    
    std::vector<Reference> references;
    std::unordered_map<std::string, int> referenceIndecesById;
    std::vector<std::vector<PositionHash>> positionHashesByReference;
    std::unordered_map<hash_t, std::vector<Locus>> lociByHash;
    
    Parameters parameters;
    double kmerSpace;
    std::string file;
};

void addMinHashes(HashSet & lociByHash, HashPriorityQueue & minHashesQueue, bloom_filter * bloomFilter, char * seq, uint32_t length, const Sketch::Parameters & parameters, uint64_t & kmersTotal, uint64_t & kmersUsed);
void getMinHashPositions(std::vector<Sketch::PositionHash> & loci, char * seq, uint32_t length, const Sketch::Parameters & parameters, int verbosity = 0);
bool hasSuffix(std::string const & whole, std::string const & suffix);
void reverseComplement(const char * src, char * dest, int length);

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);

#endif
