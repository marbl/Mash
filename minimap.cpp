#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "MurmurHash3.h"
#include <map>
#include <vector>
#include <unordered_map>

using namespace::std;

KSEQ_INIT(gzFile, gzread)

const static int kmer = 4;
const static int seed = 42; // TODO: better seed???
const static int mins = 5;

typedef uint32_t hash_t;

struct Locus
{
	Locus(uint32_t sequenceNew, uint32_t positionNew) :
		sequence(sequenceNew),
		position(positionNew)
		{}
	
	uint32_t sequence;
	uint32_t position;
};

typedef map < hash_t, vector<Locus> > LociByHash_map;
typedef unordered_map < hash_t, vector<Locus> > LociByHash_umap;

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
	
	LociByHash_umap lociByHashGlobal;
	
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
		LociByHash_map lociByHashLocal;
		
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
		
		for ( int i = 0; i < seq->seq.l - kmer + 1; i++ )
		{
			hash_t hash;
			MurmurHash3_x86_32(seq->seq.s + i, kmer, seed, &hash);
			printf("   Hash at pos %d:\t%u\n", i, hash);
			
			if
			(
				lociByHashLocal.count(hash) == 0 &&
				(
					lociByHashLocal.size() < mins ||
					hash < lociByHashLocal.rbegin()->first
				)
			)
			{
				lociByHashLocal.insert(pair<hash_t, vector<Locus> >(hash, vector<Locus>())); // insert empty vector
				
				if ( lociByHashLocal.size() > mins )
				{
					lociByHashLocal.erase(--lociByHashLocal.end());
				}
			}
			
			if ( lociByHashLocal.count(hash) )
			{
				lociByHashLocal[hash].push_back(Locus(count, i));
			}
		}
		
		printf("\n");
		
		for ( LociByHash_map::iterator i = lociByHashLocal.begin(); i != lociByHashLocal.end(); i++ )
		{
			printf("Hash %u:\n", i->first);
		
			for ( int j = 0; j < i->second.size(); j++ )
			{
				printf("   Seq: %d\tPos: %d\n", i->second.at(j).sequence, i->second.at(j).position);
			}
			
			const vector<Locus> & lociLocal = i->second;
			vector<Locus> & lociGlobal = lociByHashGlobal[i->first]; // creates if needed
			
			lociGlobal.insert(lociGlobal.end(), lociLocal.begin(), lociLocal.end());
		}
		
		count++;
	}
	
	if ( l != -1 ) printf("ERROR: return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);
	
	printf("\n");
	
	for ( LociByHash_umap::iterator i = lociByHashGlobal.begin(); i != lociByHashGlobal.end(); i++ )
	{
		printf("Hash %u:\n", i->first);
		
		for ( int j = 0; j < i->second.size(); j++ )
		{
			printf("   Seq: %d\tPos: %d\n", i->second.at(j).sequence, i->second.at(j).position);
		}
	}
	
	return 0;
}
