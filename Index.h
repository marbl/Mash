#include <capnp/message.h>
#include <capnp/serialize.h>
#include "MinHash.capnp.h"
#include <unordered_map>
#include <vector>
#include <string>

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

	int initFromCapnp(const char * file);
	int initFromSequence(const std::vector<std::string> & files, int kmerSizeNew, float compressionFactorNew);
	int writeToCapnp(const char * file) const;
	
private:
	
	std::vector<Reference> references;
	LociByHash_umap lociByHash;
	int kmerSize;
	float compressionFactor;
};

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);
