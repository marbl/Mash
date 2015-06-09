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
    useOption("sketchSize");
    useOption("concat");
    useOption("noncanonical");
    useOption("unique");
    useOption("genome");
    useOption("memory");
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
    bool list = options.at("list").active;
    
    Sketch::Parameters parameters;
    
    parameters.kmerSize = options.at("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = options.at("sketchSize").getArgumentAsNumber();
    parameters.concatenated = options.at("concat").active;
    parameters.noncanonical = options.at("noncanonical").active;
    parameters.bloomFilter = options.at("unique").active;
    parameters.genomeSize = options.at("genome").getArgumentAsNumber();
    parameters.memoryMax = options.at("memory").getArgumentAsNumber();
    
    if ( options.at("genome").active || options.at("memory").active )
    {
        parameters.bloomFilter = true;
    }
    
    if ( parameters.bloomFilter )
    {
        parameters.concatenated = true;
    }
    
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
        
        parameters.kmerSize = sketch.getKmerSize();
        parameters.noncanonical = sketch.getNoncanonical();
    }
    else
    {
        bool sketchFileExists = sketch.initHeaderFromBaseIfValid(fileReference, false);
        
        if
        (
            (options.at("kmer").active && parameters.kmerSize != sketch.getKmerSize())
        )
        {
            sketchFileExists = false;
        }
        
        if ( false && sketchFileExists )
        {
            sketch.initFromBase(fileReference, false);
            parameters.kmerSize = sketch.getKmerSize();
            parameters.noncanonical = sketch.getNoncanonical();
        }
        else
        {
            vector<string> refArgVector;
            refArgVector.push_back(fileReference);
            
            //cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;
            cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...\n";
            
            sketch.initFromSequence(refArgVector, parameters);
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
        
        threadPool.runWhenThreadAvailable(new CompareInput(sketch, sketchQuery, queryFiles[i], parameters));
        
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
        
        sketchQuery->initFromSequence(fileVector, data->parameters);
    }
    
    output->pairs.resize(sketchRef.getReferenceCount() * sketchQuery->getReferenceCount());
    
    for ( int i = 0; i < sketchQuery->getReferenceCount(); i++ )
    {
        for ( int j = 0; j < sketchRef.getReferenceCount(); j++ )
        {
            int pairIndex = i * sketchRef.getReferenceCount() + j;
            
            output->pairs[pairIndex].score = compareSketches(sketchRef.getReference(j).hashesSorted, sketchQuery->getReference(i).hashesSorted, sketchRef.getMinHashesPerWindow(), sketchQuery->getMinHashesPerWindow());
            output->pairs[pairIndex].nameRef = sketchRef.getReference(j).name;
            output->pairs[pairIndex].nameQuery = sketchQuery->getReference(i).name;
        }
    }
    
    delete data->sketchQuery;
    
    return output;
}

float compareSketches(const HashList & hashesSortedRef, const HashList & hashesSortedQuery, int targetSizeRef, int targetSizeQuery)
{
    int targetSize = targetSizeQuery < targetSizeRef ?
        targetSizeQuery :
        targetSizeRef;
    
    int i = 0;
    int j = 0;
    int common = 0;
    int denom = 0;
    
    while ( denom < targetSize && i < hashesSortedRef.size() && j < hashesSortedQuery.size() )
    {
        if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQuery.at(j), hashesSortedRef.get64()) )
        {
            i++;
        }
        else if ( hashLessThan(hashesSortedQuery.at(j), hashesSortedRef.at(i), hashesSortedRef.get64()) )
        {
            j++;
        }
        else
        {
            i++;
            j++;
            common++;
        }
        
        denom++;
    }
    
    if ( denom < targetSize )
    {
        // complete the union operation if possible
        
        if ( i < hashesSortedRef.size() )
        {
            denom += hashesSortedRef.size() - i;
        }
        
        if ( j < hashesSortedQuery.size() )
        {
            denom += hashesSortedQuery.size() - j;
        }
        
        if ( denom > targetSize )
        {
            denom = targetSize;
        }
    }
    
    return 1. - float(common) / denom;
}
