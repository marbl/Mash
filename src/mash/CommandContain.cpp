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
#include <math.h>

using namespace::std;

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
    useOption("threads");
    useOption("kmer");
    useOption("sketchSize");
    useOption("individual");
    useOption("noncanonical");
    useOption("unique");
    useOption("memory");
    useOption("genome");
    useOption("bloomError");
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
    
    parameters.kmerSize = options.at("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = options.at("sketchSize").getArgumentAsNumber();
    parameters.concatenated = ! options.at("individual").active;
    parameters.noncanonical = options.at("noncanonical").active;
    parameters.error = options.at("errorThreshold").getArgumentAsNumber();
    parameters.bloomFilter = options.at("unique").active;
    parameters.genomeSize = options.at("genome").getArgumentAsNumber();
    parameters.memoryMax = options.at("memory").getArgumentAsNumber();
    parameters.bloomError = options.at("bloomError").getArgumentAsNumber();
    
    if ( options.at("genome").active || options.at("memory").active )
    {
        parameters.bloomFilter = true;
    }
    
    if ( parameters.bloomFilter )
    {
        parameters.concatenated = true;
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
        if ( options.at("sketchSize").active )
        {
            if ( parameters.bloomFilter && parameters.minHashesPerWindow != sketchRef.getMinHashesPerWindow() )
            {
                cerr << "ERROR: The sketch size must match the reference when using a bloom filter (leave this option out to inherit from the reference sketch)." << endl;
                return 1;
            }
        }
        else
        {
            parameters.minHashesPerWindow = sketchRef.getMinHashesPerWindow();
        }
        
        parameters.kmerSize = sketchRef.getKmerSize();
        parameters.noncanonical = sketchRef.getNoncanonical();
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
    
    sketchQuery.initFromFiles(queryFiles, parameters, 0, true);
    
    for ( uint64_t i = 0; i < sketchQuery.getReferenceCount(); i++ )
    {
        for ( uint64_t j = 0; j < sketchRef.getReferenceCount(); j++ )
        {
	        threadPool.runWhenThreadAvailable(new ContainInput(sketchRef, sketchQuery, j, i, parameters));
        
			while ( threadPool.outputAvailable() )
			{
				writeOutput(threadPool.popOutputWhenAvailable(), parameters.error);
			}
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
	if ( output->error <= error )
	{
		cout << output->score << '\t' << output->error << '\t' << output->sketchRef.getReference(output->indexRef).name << '\t' << output->sketchQuery.getReference(output->indexQuery).name << endl;
	}
    
    delete output;
}

CommandContain::ContainOutput * contain(CommandContain::ContainInput * input)
{
    const Sketch & sketchRef = input->sketchRef;
    const Sketch & sketchQuery = input->sketchQuery;
    
    CommandContain::ContainOutput * output = new CommandContain::ContainOutput(input->sketchRef, input->sketchQuery, input->indexRef, input->indexQuery);
    
	output->score = containSketches(sketchRef.getReference(input->indexRef).hashesSorted, sketchQuery.getReference(input->indexQuery).hashesSorted, output->error);
    
    return output;
}

float containSketches(const HashList & hashesSortedRef, const HashList & hashesSortedQuery, float & errorToSet)
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
    
    return float(common) / j;
}
