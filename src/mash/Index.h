#ifndef Index_h
#define Index_h

#include <capnp/message.h>
#include <capnp/serialize.h>
#include "mash/capnp/MinHash.capnp.h"
#include <unordered_map>
#include <vector>
#include <string>
#include <unordered_set>

static const int seed = 42; // TODO: better seed???

static const char * capnpHeader = "Cap'n Proto";
static const int capnpHeaderLength = strlen(capnpHeader);

class Index
{
public:
    
    struct Locus
    {
        Locus(uint32_t sequenceNew, uint32_t positionNew) :
            sequence(sequenceNew),
            position(positionNew)
            {}

        uint32_t sequence;
        uint32_t position;
    };
    
    struct Reference
    {
        // no sequence for now
        
        std::string name;
        std::string comment;
        uint32_t length;
    };
    
    typedef uint32_t hash_t;
    typedef std::unordered_map < Index::hash_t, std::vector<Index::Locus> > LociByHash_umap;
	typedef std::unordered_set<Index::hash_t> Hash_set;
	
	float getCompressionFactor() const {return compressionFactor;}
	const LociByHash_umap & getLociByHash() const {return lociByHash;}
	int getKmerSize() const {return kmerSize;}
    int initFromCapnp(const char * file);
    int initFromSequence(const std::vector<std::string> & files, int kmerSizeNew, float compressionFactorNew);
    int writeToCapnp(const char * file) const;
    
private:
    
    std::vector<Reference> references;
    LociByHash_umap lociByHash;
    int kmerSize;
    float compressionFactor;
};

void getMinHashes(Index::Hash_set & lociByHash, char * seq, uint32_t length, uint32_t seqId, int kmerSize, float compressionFactor);
void getMinHashPositions(Index::LociByHash_umap & lociByHash, char * seq, uint32_t length, uint32_t seqId, int kmerSize, float compressionFactor);

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);

#endif
