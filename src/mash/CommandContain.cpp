#include "CommandContain.h"
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"

using namespace::std;

CommandContain::CommandContain()
: Command()
{
    name = "within";
    description = "Compute the containment of each query sequence (or file with -f) in the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or mash sketch files (.msh) with matching kmer sizes (-k). The score is the number of intersecting min-hashes divided by the query set size.";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    useOption("threads");
    useOption("kmer");
    useOption("sketchSize");
    useOption("concat");
    useOption("noncanonical");
    addOption("list", Option(Option::Boolean, "l", "Query files are lists of file names.", ""));
}

int CommandContain::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    int kmerSize = options.at("kmer").getArgumentAsNumber();
    int sketchSize = options.at("sketchSize").getArgumentAsNumber();
    bool concat = options.at("concat").active;
    bool list = options.at("list").active;
    bool noncanonical = options.at("noncanonical").active;
    
    Sketch sketch;
    
    const string & fileReference = arguments[0];
    
    if ( hasSuffix(fileReference, suffixSketch) )
    {
        if ( options.at("kmer").active )
        {
            cerr << "ERROR: The option " << options.at("kmer").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch.\n";
            return 1;
        }
        
        if ( options.at("noncanonical").active )
        {
            cerr << "ERROR: The option " << options.at("noncanonical").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch.\n";
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
            
            sketch.initFromSequence(refArgVector, kmerSize, sketchSize, false, 0, concat, noncanonical);
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
    
    ThreadPool<ContainInput, ContainOutput> threadPool(contain, threads);
    
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
            
            if ( sketchQuery->getNoncanonical() != sketch.getNoncanonical() )
            {
                cerr << "\nWARNING: The query sketch " << queryFiles[i] << " is " << (sketchQuery->getNoncanonical() ? "noncanonical" : "canonical") << " but the reference sketch is not. This query will be skipped.\n\n";
                delete sketchQuery;
                continue;
            }
            
            // init fully
            //
            sketchQuery->initFromCapnp(queryFiles[i].c_str());
        }
        
        threadPool.runWhenThreadAvailable(new ContainInput(sketch, sketchQuery, queryFiles[i], sketch.getKmerSize(), sketchSize, concat, sketch.getNoncanonical()));
        
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

void CommandContain::writeOutput(ContainOutput * output) const
{
    for ( int i = 0; i < output->pairs.size(); i++ )
    {
        cout << output->pairs.at(i).score << '\t' << output->pairs.at(i).nameRef << '\t' << output->pairs.at(i).nameQuery << endl;
    }
    
    delete output;
}

CommandContain::ContainOutput * contain(CommandContain::ContainInput * data)
{
    const Sketch & sketchRef = data->sketchRef;
    Sketch * sketchQuery = data->sketchQuery;
    
    CommandContain::ContainOutput * output = new CommandContain::ContainOutput();
    
    if ( sketchQuery->getReferenceCount() == 0 )
    {
        // input was sequence file; sketch now
        
        vector<string> fileVector;
        fileVector.push_back(data->file);
        
        sketchQuery->initFromSequence(fileVector, data->kmerSize, data->sketchSize, false, 0, data->concat, data->noncanonical);
    }
    
    output->pairs.resize(sketchRef.getReferenceCount() * sketchQuery->getReferenceCount());
    
    for ( int i = 0; i < sketchQuery->getReferenceCount(); i++ )
    {
        for ( int j = 0; j < sketchRef.getReferenceCount(); j++ )
        {
            int pairIndex = i * sketchRef.getReferenceCount() + j;
            
            output->pairs[pairIndex].score = containSketches(sketchRef.getReference(j).hashesSorted, sketchQuery->getReference(i).hashesSorted);
            output->pairs[pairIndex].nameRef = sketchRef.getReference(j).name;
            output->pairs[pairIndex].nameQuery = sketchQuery->getReference(i).name;
        }
    }
    
    delete data->sketchQuery;
    
    return output;
}

float containSketches(const vector<Sketch::hash_t> & hashesSortedRef, const vector<Sketch::hash_t> & hashesSortedQuery)
{
    int common = 0;
    int denom = hashesSortedRef.size() < hashesSortedQuery.size() ?
        hashesSortedRef.size() :
        hashesSortedQuery.size();
    
    int i = 0;
    int j = 0;
    
    for ( int steps = 0; steps < denom && i < hashesSortedRef.size(); steps++ )
    {
        if ( hashesSortedRef.at(i) < hashesSortedQuery.at(j) )
        {
            i++;
            steps--;
        }
        else if ( hashesSortedRef.at(i) > hashesSortedQuery.at(j) )
        {
            j++;
        }
        else
        {
            i++;
            j++;
            common++;
        }
    }
    
    return float(common) / j;
}
