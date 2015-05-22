#ifndef Sketch_h
#define Sketch_h

#include <capnp/message.h>
#include <capnp/serialize.h>
#include "mash/capnp/MinHash.capnp.h"
#include <unordered_map>
#include <vector>
#include <string>
#include <unordered_set>
#include <queue>

static const int seed = 42; // TODO: better seed???

static const char * capnpHeader = "Cap'n Proto";
static const int capnpHeaderLength = strlen(capnpHeader);

static const char * suffixSketch = ".msh";
static const char * suffixSketchWindowed = ".msw";

class Sketch
{
public:
    
#ifdef ARCH_32
    typedef uint32_t hash_t;
#else
    typedef uint64_t hash_t;
#endif
    
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
    
    typedef std::unordered_map < Sketch::hash_t, std::vector<Sketch::PositionHash> > LociByHash_umap;
    typedef std::unordered_set<Sketch::hash_t> Hash_set;
    
    struct Reference
    {
        // no sequence for now
        
        std::string name;
        std::string comment;
        uint32_t length;
        std::vector<hash_t> hashesSorted;
    };
    
    bool getConcatenated() const {return concatenated;}
    float getError() const {return error;}
    int getHashCount() const {return lociByHash.size();}
    const std::vector<Locus> & getLociByHash(hash_t hash) const;
	void getMinHashesSubsetByReference(int reference, int countToGet, Hash_set & minHashesToFill) const;
    float getMinHashesPerWindow() const {return minHashesPerWindow;}
    const Reference & getReference(int index) const {return references.at(index);}
    int getReferenceCount() const {return references.size();}
    int getReferenceIndex(std::string id) const;
    int getKmerSize() const {return kmerSize;}
    int getWindowSize() const {return windowSize;}
    bool getNoncanonical() const {return noncanonical;}
    bool hasLociByHash(hash_t hash) const {return lociByHash.count(hash);}
    void initFromBase(const std::string & file, bool windowed);
    int initFromCapnp(const char * file, bool headerOnly = false);
    int initFromSequence(const std::vector<std::string> & files, int kmerSizeNew, int sketchSizeNew, bool windowedNew, int windowSizeNew, bool concat, bool noncanonicalNew, int verbosity = 0);
    bool initHeaderFromBaseIfValid(const std::string & file, bool windowed);
    bool writeToFile() const;
    int writeToCapnp(const char * file) const;
    
private:
    
    void createIndex();
	void setMinHashesForReference(int referenceIndex, const Hash_set & hashes);
    
    std::vector<Reference> references;
    std::unordered_map<std::string, int> referenceIndecesById;
    std::vector<std::vector<PositionHash>> positionHashesByReference;
    std::unordered_map<hash_t, std::vector<Locus>> lociByHash;
    
    int kmerSize;
    float error;
    int minHashesPerWindow;
    int windowSize;
    bool windowed;
    bool concatenated;
    bool noncanonical;
    std::string file;
};

void addMinHashes(Sketch::Hash_set & lociByHash, std::priority_queue<Sketch::hash_t> & minHashesQueue, char * seq, uint32_t length, int kmerSize, int mins, bool noncanonical = false);
Sketch::hash_t getHash(const char * seq, int length);
void getMinHashPositions(std::vector<Sketch::PositionHash> & loci, char * seq, uint32_t length, int kmerSize, int minHashesPerWindow, int windowSize, int verbosity = 0);
bool hasSuffix(std::string const & whole, std::string const & suffix);
void reverseComplement(const char * src, char * dest, int length);

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);

#endif
