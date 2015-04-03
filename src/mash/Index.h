#ifndef Index_h
#define Index_h

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

static const char * suffix = ".msh";
static const char * suffixWindowed = ".msw";

class Index
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
    
    typedef std::unordered_map < Index::hash_t, std::vector<Index::PositionHash> > LociByHash_umap;
    typedef std::unordered_set<Index::hash_t> Hash_set;
    
    struct Reference
    {
        // no sequence for now
        
        std::string name;
        std::string comment;
        uint32_t length;
        Hash_set hashes;
    };
    
    bool getConcatenated() const {return concatenated;}
    int getHashCount() const {return lociByHash.size();}
    const std::vector<Locus> & getLociByHash(hash_t hash) const;
    float getMinHashesPerWindow() const {return minHashesPerWindow;}
    const Reference & getReference(int index) const {return references.at(index);}
    int getReferenceCount() const {return references.size();}
    int getReferenceIndex(std::string id) const;
    int getKmerSize() const {return kmerSize;}
    int getWindowSize() const {return windowSize;}
    bool hasLociByHash(hash_t hash) const {return lociByHash.count(hash);}
    void initFromBase(const std::string & file, bool windowed);
    int initFromCapnp(const char * file, bool headerOnly = false);
    int initFromSequence(const std::vector<std::string> & files, int kmerSizeNew, int minHashesPerWindowNew, bool windowedNew, int windowSizeNew, bool concat, int verbosity = 0);
    bool initHeaderFromBaseIfValid(const std::string & file, bool windowed);
    bool writeToFile() const;
    int writeToCapnp(const char * file) const;
    
private:
    
    void createTables();
    
    std::vector<Reference> references;
    std::unordered_map<std::string, int> referenceIndecesById;
    std::vector<std::vector<PositionHash>> positionHashesByReference;
    std::unordered_map<hash_t, std::vector<Locus>> lociByHash;
    
    int kmerSize;
    int minHashesPerWindow;
    int windowSize;
    bool windowed;
    bool concatenated;
    std::string file;
};

void addMinHashes(Index::Hash_set & lociByHash, std::priority_queue<Index::hash_t> & minHashesQueue, char * seq, uint32_t length, int kmerSize, int mins);
Index::hash_t getHash(const char * seq, int length);
void getMinHashPositions(std::vector<Index::PositionHash> & loci, char * seq, uint32_t length, int kmerSize, int minHashesPerWindow, int windowSize, int verbosity = 0);

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);

#endif
