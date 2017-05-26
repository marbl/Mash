// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandContain.h"
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"
#include "sketchParameterSetup.h"
#include <math.h>

using namespace::std;

namespace mash {

CommandContain::CommandContain()
: Command()
{
    name = "within";
    summary = "Estimate the containment of query sequences within references.";
    description = "Estimate the containment of each query file (or sequence with -i) in the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or mash sketch files (.msh) with matching k-mer sizes. Query files can also be files of file names (see -l). The score is the number of intersecting min-hashes divided by the query set size. The output format is [score, error-bound, reference-ID, query-ID].";
    argumentString = "<reference> <query> [<query>] ...";
    
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Each query file contains a list of sequence files, one per line. The reference file is not affected.", ""));
    addOption("errorThreshold", Option(Option::Number, "e", "Output", "Error bound threshold for reporting scores values. Error bounds can generally be increased by increasing the sketch size of the reference.", "0.05"));
    useOption("help");
    useSketchOptions();
}

int CommandContain::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    bool list = options.at("list").active;
    
    Sketch::Parameters parameters;
    
    parameters.error = options.at("errorThreshold").getArgumentAsNumber();
    
    if ( sketchParameterSetup(parameters, *(Command *)this) )
    {
    	return 1;
    }
    
    Sketch sketchRef;
    
    const string & fileReference = arguments[0];
    
    bool isSketch = hasSuffix(fileReference, suffixSketch);
    
    if ( isSketch )
    {
        if ( options.at("kmer").active )
        {
            cerr << "ERROR: The option " << options.at("kmer").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
            return 1;
        }
        
        if ( options.at("noncanonical").active )
        {
            cerr << "ERROR: The option " << options.at("noncanonical").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
            return 1;
        }
        
        if ( options.at("protein").active )
        {
            cerr << "ERROR: The option " << options.at("protein").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
        }
        
        if ( options.at("alphabet").active )
        {
            cerr << "ERROR: The option " << options.at("alphabet").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
        }
    }
    else
    {
        cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...";
    }
    
    vector<string> refArgVector;
    refArgVector.push_back(fileReference);
    
    //cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;
    
    sketchRef.initFromFiles(refArgVector, parameters);
    
    if ( isSketch )
    {
        parameters.minHashesPerWindow = sketchRef.getMinHashesPerWindow();
        parameters.kmerSize = sketchRef.getKmerSize();
        parameters.noncanonical = sketchRef.getNoncanonical();
        parameters.preserveCase = sketchRef.getPreserveCase();
        parameters.seed = sketchRef.getHashSeed();
        
        string alphabet;
        sketchRef.getAlphabetAsString(alphabet);
        setAlphabetFromString(parameters, alphabet.c_str());
    }
    else
    {
	    cerr << "done.\n";
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
    
    Sketch sketchQuery;
    
    sketchQuery.initFromFiles(queryFiles, parameters, 0, true, true);
    
    uint64_t pairCount = sketchRef.getReferenceCount() * sketchQuery.getReferenceCount();
    uint64_t pairsPerThread = pairCount / parameters.parallelism;
    
    if ( pairsPerThread == 0 )
    {
    	pairsPerThread = 1;
    }
    
    static uint64_t maxPairsPerThread = 0x1000;
    
    if ( pairsPerThread > maxPairsPerThread )
    {
        pairsPerThread = maxPairsPerThread;
    }
    
    uint64_t iFloor = pairsPerThread / sketchRef.getReferenceCount();
    uint64_t iMod = pairsPerThread % sketchRef.getReferenceCount();
    
    for ( uint64_t i = 0, j = 0; i < sketchQuery.getReferenceCount(); i += iFloor, j += iMod )
    {
        if ( j >= sketchRef.getReferenceCount() )
        {
            if ( i == sketchQuery.getReferenceCount() - 1 )
            {
                break;
            }
            
            i++;
            j -= sketchRef.getReferenceCount();
        }
        
		threadPool.runWhenThreadAvailable(new ContainInput(sketchRef, sketchQuery, j, i, pairsPerThread, parameters));
	
		while ( threadPool.outputAvailable() )
		{
			writeOutput(threadPool.popOutputWhenAvailable(), parameters.error);
		}
    }
    
    while ( threadPool.running() )
    {
        writeOutput(threadPool.popOutputWhenAvailable(), parameters.error);
    }
    
    return 0;
}

void CommandContain::writeOutput(ContainOutput * output, float error) const
{
    uint64_t i = output->indexQuery;
    uint64_t j = output->indexRef;
    
    for ( uint64_t k = 0; k < output->pairCount && i < output->sketchQuery.getReferenceCount(); k++ )
    {
        const ContainOutput::PairOutput * pair = &output->pairs[k];
        
		if ( pair->error <= error )
		{
			cout << pair->score << '\t' << pair->error << '\t' << output->sketchRef.getReference(j).name << '\t' << output->sketchQuery.getReference(i).name << endl;
		}
        
        j++;
        
        if ( j == output->sketchRef.getReferenceCount() )
        {
            j = 0;
            i++;
        }
	}
    
    delete output;
}

CommandContain::ContainOutput * contain(CommandContain::ContainInput * input)
{
    const Sketch & sketchRef = input->sketchRef;
    const Sketch & sketchQuery = input->sketchQuery;
    
    CommandContain::ContainOutput * output = new CommandContain::ContainOutput(input->sketchRef, input->sketchQuery, input->indexRef, input->indexQuery, input->pairCount);
    
    uint64_t i = input->indexQuery;
    uint64_t j = input->indexRef;
    
    for ( uint64_t k = 0; k < input->pairCount && i < sketchQuery.getReferenceCount(); k++ )
    {
		output->pairs[k].score = containSketches(sketchRef.getReference(j).hashesSorted, sketchQuery.getReference(i).hashesSorted, output->pairs[k].error);
        
        j++;
        
        if ( j == sketchRef.getReferenceCount() )
        {
            j = 0;
            i++;
        }
    }
    
    return output;
}

double containSketches(const HashList & hashesSortedRef, const HashList & hashesSortedQuery, double & errorToSet)
{
    int common = 0;
    int denom = hashesSortedRef.size() < hashesSortedQuery.size() ?
        hashesSortedRef.size() :
        hashesSortedQuery.size();
    
    int i = 0;
    int j = 0;
    
    for ( int steps = 0; steps < denom && i < hashesSortedRef.size(); steps++ )
    {
        if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQuery.at(j), hashesSortedRef.get64()) )
        {
            i++;
            steps--;
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
    }
    
    errorToSet = 1. / sqrt(j);
    
    return double(common) / j;
}

} // namespace mash
