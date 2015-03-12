#include "CommandFind.h"
#include "Index.h"
#include <zlib.h>
#include "kseq.h"
#include <iostream>
#include <set>
#include "ThreadPool.h"

using namespace::std;

KSEQ_INIT(gzFile, gzread)

CommandFind::CommandFind()
{
    name = "find";
    description = "Compare query sequences to a reference index";
    argumentString = "index.mash fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
    addOption("threshold", Option(Option::Number, "t", "Threshold. This fraction of the query sequence's min-hashes must appear in a query-sized window of a reference sequence for the match to be reported.", "0.2"));
    addOption("threads", Option(Option::Number, "p", "Parallelism. This many threads will be spawned to perform the find, each one handling on query sequence at a time.", "1"));
}

int CommandFind::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    float threshold = options.at("threshold").getArgumentAsNumber(0, 1);
    int threads = options.at("threads").getArgumentAsNumber();
    
    Index index;
    //
    index.initFromCapnp(arguments[0].c_str());
    
    int l;
    int count = 0;
    
    ThreadPool<FindInput, FindOutput> threadPool(find, threads);
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        gzFile fp = gzopen(arguments[i].c_str(), "r");
        kseq_t *seq = kseq_init(fp);
        
        while ((l = kseq_read(seq)) >= 0)
        {
            if ( l < index.getKmerSize() )
            {
                continue;
            }
            
            //printf("Query name: %s\tlength: %d\n\n", seq->name.s, l);
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            
            threadPool.runWhenThreadAvailable(new FindInput(index, seq->name.s, seq->seq.s, l, threshold));
            
            while ( threadPool.outputAvailable() )
            {
                //cout << "popping\n";
                writeOutput(index, threadPool.popOutputWhenAvailable());
            }
        }
        
        if ( l != -1 )
        {
            printf("ERROR: return value: %d\n", l);
            return 1;
        }
        
        kseq_destroy(seq);
        gzclose(fp);
    }
    
    while ( threadPool.running() )
    {
        writeOutput(index, threadPool.popOutputWhenAvailable());
    }
    
    return 0;
}

void CommandFind::writeOutput(const Index & index, const FindOutput * output) const
{
    //cout << output->seqId << endl;
    for ( int i = 0; i < output->hits.size(); i++ )
    {
        const FindOutput::Hit & hit = output->hits.at(i);
        
        //cout << output->seqId << '\t' << index.getReference(hit.ref).name << '\t' << hit.start << '\t' << hit.end << '\t' << hit.score << endl;
    }
    
    delete output;
}

CommandFind::FindOutput * find(CommandFind::FindInput * data)
{
    CommandFind::FindOutput * output = new CommandFind::FindOutput();
    
    typedef unordered_map < uint32_t, set<uint32_t> > PositionsBySequence_umap;
    
    Index::Hash_set minHashes;
    
    const Index & index = data->index;
    int kmerSize = index.getKmerSize();
    float compressionFactor = index.getCompressionFactor();
    int windowSize = index.getWindowSize();
    int length = data->length;
    char * seq = data->seq;
    float threshold = data->threshold;
    
    output->seqId = data->seqId;
    
    int mins = length / compressionFactor;
    //
    if ( mins < 1 )
    {
        mins = 1;
    }
    
    //cout << "Mins: " << mins << "\t length: " << length << "\tComp: " << compressionFactor << endl;
    getMinHashes(minHashes, seq, length, 0, kmerSize, compressionFactor);
    
    for ( int i = 0; i < index.getLociByReference().size(); i++ )
    {
        const vector<Index::Locus> & loci = index.getLociByReference().at(i);
        
        // pointer to the position at the beginning of the window; to be updated
        // as the end of the window is incremented
        //
        vector<Index::Locus>::const_iterator windowStart = loci.begin();
        
        //cout << "Clustering in seq " << i->first << endl;
        
        // the number of positions between the window start and end (inclusive)
        //
        int windowCount = 0;
        
        for ( vector<Index::Locus>::const_iterator j = loci.begin(); j != loci.end(); j++ )
        {
            windowCount++;
            
            //cout << *windowStart << "\t" << *j << endl;
            
            // update window start if it is too far behind
            //
            while ( windowStart != j && j->position > length && windowStart->position < j->position - length + 1 )
            {
                //cout << "moving " << *j - length + 1 << endl;
                windowStart++;
                windowCount--;
            }
            
            // extend the right of the window if possible
            //
            while ( j != loci.end() && j->position - windowStart->position < length )
            {
                windowCount++;
                j++;
            }
            //
            windowCount--;
            j--;
            
            //cout << *windowStart << "\t" << *j << endl;
            float score = float(windowCount) / mins;
            
            if ( score >= threshold )
            {
                cout << data->seqId << '\t' << index.getReference(i).name << '\t' << windowStart->position << '\t' << j->position << '\t' << float(windowCount) / mins << endl;
                
                output->hits.resize(output->hits.size() + 1);
                
                CommandFind::FindOutput::Hit & hit = output->hits[output->hits.size() - 1];
                
                hit.ref = i;
                hit.start = windowStart->position;
                hit.end = j->position;
                hit.score = score;
                
                for ( vector<Index::Locus>::const_iterator k = windowStart; k != loci.end() && k->position <= j->position; k++ )
                {
                    cout << "      " << k->position << endl;
                }
            }
        }
    }
    //cout << "done\n";
    return output;
}
