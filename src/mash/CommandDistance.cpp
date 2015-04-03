#include "CommandDistance.h"
#include "Index.h"
#include <iostream>
#include <zlib.h>
#include "kseq.h"
#include "ThreadPool.h"

using namespace::std;

KSEQ_INIT(gzFile, gzread)

CommandDistance::CommandDistance()
: Command()
{
    name = "dist";
    description = "Compute the global distance of each input sequence to the reference index. The score is computed as the Jaccard index, which is the intesection divided by the union, for the sets of min-hashes in the reference and query.";
    argumentString = "reference.mash fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
    addOption("threads", Option(Option::Number, "p", "Parallelism. This many threads will be spawned, each one handling one input file at a time.", "1"));
    addOption("kmer", Option(Option::Number, "k", "Kmer size. Hashes will be based on strings of this many nucleotides.", "11"));
    addOption("mins", Option(Option::Number, "m", "Min-hashes per input file.", "10000"));
    addOption("concat", Option(Option::Boolean, "c", "Concatenate multi-fasta seqeunces (and use files for names). Kmers across boundaries will not be considered. Not compatible with -w.", ""));
}

int CommandDistance::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    int kmerSize = options.at("kmer").getArgumentAsNumber();
    int mins = options.at("mins").getArgumentAsNumber();
    bool concat = options.at("concat").active;
    
    Index index;
    
    bool indexFileExists = index.initHeaderFromBaseIfValid(arguments[0], false);
    
    if
    (
        (options.at("kmer").active && kmerSize != index.getKmerSize()) ||
        (options.at("mins").active && mins != index.getMinHashesPerWindow()) ||
        concat != index.getConcatenated()
    )
    {
        indexFileExists = false;
    }
    
    if ( indexFileExists )
    {
        index.initFromBase(arguments[0], false);
    }
    else
    {
        vector<string> refArgVector;
        refArgVector.push_back(arguments[0]);
        
        cerr << "Sketch for " << arguments[0] << " not found or out of date; creating..." << endl;
        index.initFromSequence(refArgVector, kmerSize, mins, false, 0, concat);
        
        if ( index.writeToFile() )
        {
            cerr << "Sketch saved for subsequent runs." << endl;
        }
        else
        {
            cerr << "The sketch for " << arguments[0] << " could not be saved; it will be sketched again next time." << endl;
        }
    }
    
    ThreadPool<CompareInput, CompareOutput> threadPool(compare, threads);
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        for ( int j = 0; j < index.getReferenceCount(); j++ )
        {
            threadPool.runWhenThreadAvailable(new CompareInput(index.getReference(j).hashes, index.getReference(j).name, arguments[i], kmerSize, mins, concat));
        
            while ( threadPool.outputAvailable() )
            {
                writeOutput(threadPool.popOutputWhenAvailable());
            }
        }
    }
    
    while ( threadPool.running() )
    {
        writeOutput(threadPool.popOutputWhenAvailable());
    }
    
    return 0;
}

void CommandDistance::writeOutput(CompareOutput * output) const
{
    for ( int i = 0; i < output->pairs.size(); i++ )
    {
        cout << output->pairs.at(i).score << '\t' << output->nameRef << '\t' << output->pairs.at(i).file << endl;
    }
    
    delete output;
}

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * data)
{
    const Index::Hash_set & minHashesRef = data->minHashesRef;
    const string file = data->file;
    
    CommandDistance::CompareOutput * output = new CommandDistance::CompareOutput();
    Index index;
    vector<string> fileVector;
    fileVector.push_back(file);
    
    index.initFromSequence(fileVector, data->kmerSize, data->mins, false, 0, data->concat);
    output->nameRef = data->nameRef;
    output->pairs.resize(index.getReferenceCount());
    
    for ( int i = 0; i < index.getReferenceCount(); i++ )
    {
        int common = 0;
        
        const Index::Hash_set & minHashes = index.getReference(i).hashes;
        
        for ( Index::Hash_set::const_iterator i = minHashes.begin(); i != minHashes.end(); i++ )
        {
            if ( minHashesRef.count(*i) == 1 )
            {
                common++;
            }
        }
        
		int denominator;
		
		if ( minHashesRef.size() >= data->mins & minHashes.size() >= data->mins )
		{
			denominator = data->mins;
		}
		else
		{
			if ( minHashesRef.size() > minHashes.size() )
			{
				denominator = minHashesRef.size();
			}
			else
			{
				denominator = minHashes.size();
			}
		}
		
        output->pairs[i].score = float(common) / denominator;
        output->pairs[i].file = index.getReference(i).name;
    }
    
    return output;
}
