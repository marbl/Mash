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
    addOption("kmer", Option(Option::Number, "k", "Kmer size. Hashes will be based on strings of this many nucleotides.", "11"));
    addOption("mins", Option(Option::Number, "m", "Min-hashes per input file.", "10000"));
}

int CommandCompare::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    int kmerSize = options.at("kmer").getArgumentAsNumber();
    int mins = options.at("mins").getArgumentAsNumber();
    
    Index::Hash_set minHashesRef;
    getMinHashesForFile(minHashesRef, arguments[0], kmerSize, mins);
    
    ThreadPool<CompareInput, CompareOutput> threadPool(compare, threads);
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        threadPool.runWhenThreadAvailable(new CompareInput(minHashesRef, arguments[i], kmerSize, mins));
        
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
    const Index::Hash_set & minHashesRef = data->minHashesRef;
    const string file = data->file;
    
    CommandCompare::CompareOutput * output = new CommandCompare::CompareOutput();
    int common = 0;
    Index::Hash_set minHashes;
    
    getMinHashesForFile(minHashes, file, data->kmerSize, data->mins);
    
    for ( Index::Hash_set::const_iterator i = minHashes.begin(); i != minHashes.end(); i++ )
    {
        if ( minHashesRef.count(*i) == 1 )
        {
            common++;
        }
    }
    
    output->score = float(common) / (minHashesRef.size() + minHashes.size() - common);
    output->file = file;
    
    return output;
}

void getMinHashesForFile(Index::Hash_set & minHashes, const string & file, int kmerSize, int mins)
{
    int l;
    
    gzFile fp = gzopen(file.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    
    minHashes.clear();
    priority_queue<Index::hash_t> minHashesQueue;
    
    while ((l = kseq_read(seq)) >= 0)
    {
        if ( l < kmerSize )
        {
            continue;
        }
        
        //printf("Query name: %s\tlength: %d\n\n", seq->name.s, l);
        //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
        //printf("seq: %s\n", seq->seq.s);
        //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        
        //getMinHashes(minHashesLocal, seq->seq.s, l, 0, indexRef.getKmerSize(), indexRef.getMinHashesPerWindow());
        //
        addMinHashes(minHashes, minHashesQueue, seq->seq.s, l, kmerSize, mins);
    }
    
    if ( l != -1 )
    {
        printf("ERROR: return value: %d\n", l);
        exit(1);
    }
    
    kseq_destroy(seq);
    gzclose(fp);
}
