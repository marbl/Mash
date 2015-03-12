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
    
    typedef uint32_t hash_t;
    
    struct Locus
    {
        Locus(uint32_t positionNew, uint32_t hashNew) :
            position(positionNew),
            hash(hashNew)
            {}

        uint32_t position;
        hash_t hash;
    };
    
    struct Reference
    {
        // no sequence for now
        
        std::string name;
        std::string comment;
        uint32_t length;
    };
    
    typedef std::unordered_map < Index::hash_t, std::vector<Index::Locus> > LociByHash_umap;
	typedef std::unordered_set<Index::hash_t> Hash_set;
	
	float getCompressionFactor() const {return compressionFactor;}
	const Reference & getReference(int index) const {return references.at(index);}
	const std::vector<std::vector<Locus>> & getLociByReference() const {return lociByReference;}
	int getKmerSize() const {return kmerSize;}
	int getWindowSize() const {return windowSize;}
    int initFromCapnp(const char * file);
    int initFromSequence(const std::vector<std::string> & files, int kmerSizeNew, float compressionFactorNew, int windowSizeNew, bool verbose = false);
    int writeToCapnp(const char * file) const;
    
private:
    
    std::vector<Reference> references;
    std::vector<std::vector<Locus>> lociByReference;
    
    int kmerSize;
    float compressionFactor;
    int windowSize;
};

void getMinHashes(Index::Hash_set & lociByHash, char * seq, uint32_t length, uint32_t seqId, int kmerSize, float compressionFactor);
void getMinHashPositions(std::vector<Index::Locus> & loci, char * seq, uint32_t length, int kmerSize, float compressionFactor, int windowSize, bool verbose = false);

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);

#endif
