#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "MurmurHash3.h"
#include <map>
#include <vector>
#include <unordered_map>
#include <iostream>

#include "CommandIndex.h"

using namespace::std;

KSEQ_INIT(gzFile, gzread)

typedef uint32_t hash_t;
static const int seed = 42; // TODO: better seed???
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

int CommandIndex::run() const
{
	cout << "kmer size: " << options.at("kmer").argument << endl;
	
	gzFile fp;
	kseq_t *seq;
	int l;
	
	LociByHash_umap lociByHashGlobal;
	
	fp = gzopen(arguments[0].c_str(), "r");
	seq = kseq_init(fp);
	int count = 0;
	
	int kmer = atoi(options.at("kmer").argument.c_str());
	int factor = atoi(options.at("factor").argument.c_str());
	int mins = atoi(options.at("hashes").argument.c_str());
	
	while ((l = kseq_read(seq)) >= 0)
	{
		LociByHash_map lociByHashLocal;
		
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
		
		int stride = factor; // TODO: ???
		
		for ( int i = 0; i < seq->seq.l - kmer + 1; i += stride )
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
	
	printf("\nCombined hash table:\n\n");
	
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
