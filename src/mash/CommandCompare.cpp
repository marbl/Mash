#include "CommandCompare.h"
#include "Index.h"
#include <iostream>
#include <zlib.h>
#include "kseq.h"
#include "ThreadPool.h"

using namespace::std;

KSEQ_INIT(gzFile, gzread)

CommandCompare::CommandCompare()
{
    name = "compare";
    description = "Compute the global distance of each input sequence to the reference index. The score is computed as the Jaccard index, which is the intesection divided by the union, for the sets of min-hashes in the reference and query.";
    argumentString = "reference.mash fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
    addOption("threads", Option(Option::Number, "p", "Parallelism. This many threads will be spawned, each one handling one input file at a time.", "1"));
}

int CommandCompare::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    
    Index indexRef;
    indexRef.initFromCapnp(arguments[0].c_str());
    
    ThreadPool<CompareInput, CompareOutput> threadPool(compare, threads);
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        threadPool.runWhenThreadAvailable(new CompareInput(indexRef, arguments[i]));
        
        while ( threadPool.outputAvailable() )
        {
            writeOutput(threadPool.popOutputWhenAvailable());
        }
    }
    
    while ( threadPool.running() )
    {
        writeOutput(threadPool.popOutputWhenAvailable());
    }
    
    return 0;
}

void CommandCompare::writeOutput(CompareOutput * output) const
{
    cout << output->score << '\t' << output->file << endl;
    delete output;
}

CommandCompare::CompareOutput * compare(CommandCompare::CompareInput * data)
{
    const Index & indexRef = data->indexRef;
    const string file = data->file;
    
    CommandCompare::CompareOutput * output = new CommandCompare::CompareOutput();
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
        return 0;
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
    
    output->score = float(common) / (indexRef.getLociByHash().size() + minHashesGlobal.size() - common);
    output->file = file;
    
    return output;
}
