#include "CommandCompare.h"
#include "Index.h"
#include <iostream>
#include <zlib.h>
#include "kseq.h"

using namespace::std;

KSEQ_INIT(gzFile, gzread)

CommandCompare::CommandCompare()
{
    name = "compare";
    description = "Compute the global distance of each input sequence to the reference index. The score is computed as the Jaccard index, which is the intesection divided by the union, for the sets of min-hashes in the reference and query.";
    argumentString = "reference.mash fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
}

int CommandCompare::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    Index indexRef;
    indexRef.initFromCapnp(arguments[0].c_str());
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        cout << compare(indexRef, arguments[i]) << '\t' << arguments[i] << endl;
    }
    
    return 0;
}

float compare(const Index & indexRef, const string file)
{
    int common = 0;
    int l;
    
    gzFile fp = gzopen(file.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    
    Index::Hash_set minHashesGlobal;
    
    while ((l = kseq_read(seq)) >= 0)
    {
        if ( l < indexRef.getKmerSize() )
        {
            continue;
        }
        
        //printf("Query name: %s\tlength: %d\n\n", seq->name.s, l);
        //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
        //printf("seq: %s\n", seq->seq.s);
        //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        
        Index::Hash_set minHashesLocal;
        
        getMinHashes(minHashesLocal, seq->seq.s, l, 0, indexRef.getKmerSize(), indexRef.getCompressionFactor());
        
        for ( Index::Hash_set::const_iterator i = minHashesLocal.begin(); i != minHashesLocal.end(); i++ )
        {
            minHashesGlobal.insert(*i);
        }
    }
    
    if ( l != -1 )
    {
        printf("ERROR: return value: %d\n", l);
        return 1;
    }
    
    kseq_destroy(seq);
    gzclose(fp);
    
    for ( Index::Hash_set::const_iterator i = minHashesGlobal.begin(); i != minHashesGlobal.end(); i++ )
    {
        if ( indexRef.getLociByHash().count(*i) != 0 )
        {
            common++;
        }
    }
    
    return float(common) / (indexRef.getLociByHash().size() + minHashesGlobal.size() - common);
}
