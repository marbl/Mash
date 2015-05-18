#include "CommandDistance.h"
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"

using namespace::std;

CommandDistance::CommandDistance()
: Command()
{
    name = "dist";
    description = "Compute the global distance of each query sequence to the reference. Both the reference and queries can be fasta or fastq, gzipped or not. The reference can also be a mash sketch file (.msh). If the reference is fasta/fastq and a sketch file does not exist (or is out of date), one will be written with the current options and used for future runs if possible. The score is computed as the number of matching min-hashes divided by -m, or the larger number of kmers if both sequences have fewer than -m.";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    useOption("threads");
    useOption("kmer");
    useOption("factor");
    useOption("concat");
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
    float factor = options.at("factor").getArgumentAsNumber();
    bool concat = options.at("concat").active;
    
    Sketch sketch;
    
    const string & fileReference = arguments[0];
    
    if ( hasSuffix(fileReference, suffixSketch) )
    {
        if ( options.at("kmer").active )
        {
            cerr << "ERROR: The option " << options.at("kmer").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch.\n";
            return 1;
        }
        
        sketch.initFromCapnp(fileReference.c_str());
    }
    else
    {
        bool sketchFileExists = sketch.initHeaderFromBaseIfValid(fileReference, false);
        
        if
        (
            (options.at("kmer").active && kmerSize != sketch.getKmerSize())
        )
        {
            sketchFileExists = false;
        }
        
        if ( false && sketchFileExists )
        {
            sketch.initFromBase(fileReference, false);
            kmerSize = sketch.getKmerSize();
        }
        else
        {
            vector<string> refArgVector;
            refArgVector.push_back(fileReference);
            
            //cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;
            cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...\n";
            
            sketch.initFromSequence(refArgVector, kmerSize, factor, false, 0, concat);
            /*
            if ( sketch.writeToFile() )
            {
                cerr << "Sketch saved for subsequent runs." << endl;
            }
            else
            {
                cerr << "The sketch for " << fileReference << " could not be saved; it will be sketched again next time." << endl;
            }*/
        }
    }
    
    ThreadPool<CompareInput, CompareOutput> threadPool(compare, threads);
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        // If the input is a sketch file, load in the main thread; otherwise,
        // leave it to the child. Either way, the child will delete.
        //
        Sketch * sketchQuery = new Sketch();
        
        if ( hasSuffix(arguments[i], suffixSketch) )
        {
            // init header to check params
            //
            sketchQuery->initFromCapnp(arguments[i].c_str(), true);
            
            if ( sketchQuery->getKmerSize() != sketch.getKmerSize() )
            {
                cerr << "\nWARNING: The query sketch " << arguments[i] << " has a kmer size (" << sketchQuery->getKmerSize() << ") that does not match the reference sketch (" << sketch.getKmerSize() << "). This query will be skipped.\n\n";
                continue;
            }
            
            // init fully
            //
            sketchQuery->initFromCapnp(arguments[i].c_str());
        }
        
        threadPool.runWhenThreadAvailable(new CompareInput(sketch, sketchQuery, arguments[i], sketch.getKmerSize(), sketch.getFactor(), concat));
        
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

void CommandDistance::writeOutput(CompareOutput * output) const
{
    for ( int i = 0; i < output->pairs.size(); i++ )
    {
        cout << output->pairs.at(i).score << '\t' << output->pairs.at(i).nameRef << '\t' << output->pairs.at(i).nameQuery << endl;
    }
    
    delete output;
}

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * data)
{
    const Sketch & sketchRef = data->sketchRef;
    Sketch * sketchQuery = data->sketchQuery;
    
    CommandDistance::CompareOutput * output = new CommandDistance::CompareOutput();
    
    if ( sketchQuery->getReferenceCount() == 0 )
    {
        // input was sequence file; sketch now
        
        vector<string> fileVector;
        fileVector.push_back(data->file);
        
        sketchQuery->initFromSequence(fileVector, data->kmerSize, data->factor, false, 0, data->concat);
    }
    
    output->pairs.resize(sketchRef.getReferenceCount() * sketchQuery->getReferenceCount());
    
    for ( int i = 0; i < sketchQuery->getReferenceCount(); i++ )
    {
        for ( int j = 0; j < sketchRef.getReferenceCount(); j++ )
        {
            // Keep a reference the smaller of the two hash vectors (whether ref or query)
            // and get a set of the same size from the other so they can be queried.
            
            bool refSmaller = sketchRef.getReference(j).hashesSorted.size() < sketchQuery->getReference(i).hashesSorted.size();
            
            const vector<Sketch::hash_t> & hashesSorted = refSmaller ?
                sketchRef.getReference(j).hashesSorted :
                sketchQuery->getReference(i).hashesSorted;
            
            int minHashCount = hashesSorted.size();
            Sketch::Hash_set minHashesSubset;
            
            if ( refSmaller )
            {
                sketchRef.getMinHashesSubsetByReference(j, minHashCount, minHashesSubset);
            }
            else
            {
                sketchQuery->getMinHashesSubsetByReference(i, minHashCount, minHashesSubset);
            }
            
            int common = 0;
            
            for ( int k = 0; k < minHashCount; k++ )
            {
                if ( minHashesSubset.count(hashesSorted.at(k)) == 1 )
                {
                    common++;
                }
            }
            
            int pairIndex = i * sketchRef.getReferenceCount() + j;
            
            output->pairs[pairIndex].score = float(common) / minHashCount;
            output->pairs[pairIndex].nameRef = sketchRef.getReference(j).name;
            output->pairs[pairIndex].nameQuery = sketchQuery->getReference(i).name;
        }
    }
    
    delete data->sketchQuery;
    
    return output;
}
