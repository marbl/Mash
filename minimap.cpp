#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "MurmurHash3.h"
#include <map>
#include <queue>

using namespace::std;

KSEQ_INIT(gzFile, gzread)

const static int kmer = 4;
const static int seed = 42; // TODO: better seed???
const static int mins = 5;

typedef uint32_t Hash;

struct Locus
{
	Locus(uint32_t sequenceNew, uint32_t positionNew) :
		sequence(sequenceNew),
		position(positionNew)
		{}
	
	uint32_t sequence;
	uint32_t position;
};

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
	
	priority_queue<Hash> minHashes;
	map < Hash, vector<Locus> > lociByHash;
	
	if (argc == 1)
	{
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	int count = 0;
	
	while ((l = kseq_read(seq)) >= 0)
	{
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
		
		for ( int i = 0; i < seq->seq.l - kmer + 1; i++ )
		{
			Hash hash;
			MurmurHash3_x86_32(seq->seq.s + i, kmer, seed, &hash);
			printf("   Hash at pos %d:\t%u\n", i, hash);
			
			if ( lociByHash.count(hash) == 0 && ( minHashes.size() < mins || hash < minHashes.top() ) )
			{
				minHashes.push(hash);
				lociByHash[hash]; // insert empty vector
				//printf("PUSH %u, %d\n", hash, lociByHash.count(hash));
			
				if ( minHashes.size() > mins )
				{
					//printf("POP  %u\n", minHashes.top());
					lociByHash.erase(minHashes.top());
					minHashes.pop();
				}
			}
			
			if ( lociByHash.count(hash) )
			{
				//printf("adding to %u...\n", hash);
				lociByHash[hash].push_back(Locus(count, i));
			}
		}
		
		printf("\n");
		count++;
	}
	
	if ( l != -1 ) printf("ERROR: return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);
	
	printf("\n");
	
	for ( map< Hash, vector<Locus> >::iterator i = lociByHash.begin(); i != lociByHash.end(); i++ )
	{
		printf("Hash %u:\n", i->first);
		
		for ( int j = 0; j < i->second.size(); j++ )
		{
			printf("   Seq: %d\tPos: %d\n", i->second.at(j).sequence, i->second.at(j).position);
		}
	}
	
	return 0;
}
