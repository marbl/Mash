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
    description = "Compute the distance of each query sequence (or file with -f) to the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or mash sketch files (.msh) with matching kmer sizes (-k). The distance is one minus the Jaccard score for the set of min-hashes whose size is that of the smaller sketch.";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    useOption("threads");
    useOption("kmer");
    useOption("error");
    useOption("concat");
    addOption("list", Option(Option::Boolean, "l", "Query files are lists of file names.", ""));
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
    float error = options.at("error").getArgumentAsNumber();
    bool concat = options.at("concat").active;
    bool list = options.at("list").active;
    
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
            
            sketch.initFromSequence(refArgVector, kmerSize, error, false, 0, concat);
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
    
    vector<string> queryFiles;
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        if ( list )
        {
            splitFile(arguments[i], queryFiles);
        }
        else
        {
            queryFiles.push_back(arguments[i]);
        }
    }
    
    for ( int i = 0; i < queryFiles.size(); i++ )
    {
        // If the input is a sketch file, load in the main thread; otherwise,
        // leave it to the child. Either way, the child will delete.
        //
        Sketch * sketchQuery = new Sketch();
        
        if ( hasSuffix(queryFiles[i], suffixSketch) )
        {
            // init header to check params
            //
            sketchQuery->initFromCapnp(queryFiles[i].c_str(), true);
            
            if ( sketchQuery->getKmerSize() != sketch.getKmerSize() )
            {
                cerr << "\nWARNING: The query sketch " << queryFiles[i] << " has a kmer size (" << sketchQuery->getKmerSize() << ") that does not match the reference sketch (" << sketch.getKmerSize() << "). This query will be skipped.\n\n";
                delete sketchQuery;
                continue;
            }
            
            // init fully
            //
            sketchQuery->initFromCapnp(queryFiles[i].c_str());
        }
        
        threadPool.runWhenThreadAvailable(new CompareInput(sketch, sketchQuery, queryFiles[i], sketch.getKmerSize(), sketch.getError(), concat));
        
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
        
        sketchQuery->initFromSequence(fileVector, data->kmerSize, data->error, false, 0, data->concat);
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
                sketchQuery->getMinHashesSubsetByReference(i, minHashCount, minHashesSubset);
            }
            else
            {
                sketchRef.getMinHashesSubsetByReference(j, minHashCount, minHashesSubset);
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
            
            output->pairs[pairIndex].score = 1. - float(common) / minHashCount;
            output->pairs[pairIndex].nameRef = sketchRef.getReference(j).name;
            output->pairs[pairIndex].nameQuery = sketchQuery->getReference(i).name;
        }
    }
    
    delete data->sketchQuery;
    
    return output;
}
