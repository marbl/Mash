// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandDistance.h"
#include "CommandTriangle.h"
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"
#include "sketchParameterSetup.h"
#include <math.h>

#ifdef USE_BOOST
    #include <boost/math/distributions/binomial.hpp>
    using namespace::boost::math;
#else
    #include <gsl/gsl_cdf.h>
#endif

using namespace::std;

namespace mash {

CommandTriangle::CommandTriangle()
: Command()
{
    name = "triangle";
    summary = "Estimate a lower-triangular distance matrix.";
    description = "Estimate the distance of each input sequence to every other input sequence. Outputs a lower-triangular distance matrix in relaxed Phylip format. The input sequences can be fasta or fastq, gzipped or not, or Mash sketch files (.msh) with matching k-mer sizes. Input files can also be files of file names (see -l). If more than one input file is provided, whole files are compared by default (see -i).";
    argumentString = "<seq1> [<seq2>] ...";
    
    useOption("help");
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <query> specify paths to sequence files, one per line. The reference file is not affected.", ""));
    addOption("comment", Option(Option::Boolean, "C", "Output", "Use comment fields for sequence names instead of IDs.", ""));
    addOption("edge", Option(Option::Boolean, "E", "Output", "Output edge list with fields [seq1, seq2, dist, p-val, shared-hashes].", ""));
    //addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances and divide by k-mer size to provide a better analog to phylogenetic distance. The special case of zero shared min-hashes will result in a distance of 1.", ""));
    useSketchOptions();
}

int CommandTriangle::run() const
{
    if ( arguments.size() < 1 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    bool list = options.at("list").active;
    //bool log = options.at("log").active;
    double pValueMax = 0;
    bool comment = options.at("comment").active;
    bool edge = options.at("edge").active;
    
    Sketch::Parameters parameters;
    
    if ( sketchParameterSetup(parameters, *(Command *)this) )
    {
    	return 1;
    }
    
    if ( arguments.size() == 1 )
    {
    	parameters.concatenated = false;
    }
    
    Sketch sketch;
    
    uint64_t lengthMax;
    double randomChance;
    int kMin;
    string lengthMaxName;
    int warningCount = 0;
    
    vector<string> queryFiles;
    
    for ( int i = 0; i < arguments.size(); i++ )
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
    
    sketch.initFromFiles(queryFiles, parameters);
    
    double lengthThreshold = (parameters.warning * sketch.getKmerSpace()) / (1. - parameters.warning);
    
	for ( uint64_t i = 0; i < sketch.getReferenceCount(); i++ )
	{
		uint64_t length = sketch.getReference(i).length;
		
		if ( length > lengthThreshold )
		{
			if ( warningCount == 0 || length > lengthMax )
			{
				lengthMax = length;
				lengthMaxName = sketch.getReference(i).name;
				randomChance = sketch.getRandomKmerChance(i);
				kMin = sketch.getMinKmerSize(i);
			}
		
			warningCount++;
		}
	}
    
    if ( !edge )
    {
        cout << '\t' << sketch.getReferenceCount() << endl;
        cout << (comment ? sketch.getReference(0).comment : sketch.getReference(0).name) << endl;
    }
    
    ThreadPool<TriangleInput, TriangleOutput> threadPool(compare, threads);
    
    for ( uint64_t i = 1; i < sketch.getReferenceCount(); i++ )
    {
        threadPool.runWhenThreadAvailable(new TriangleInput(sketch, i, parameters));
        
        while ( threadPool.outputAvailable() )
        {
            writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValueMax);
        }
    }
    
    while ( threadPool.running() )
    {
        writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValueMax);
    }
    
    if ( !edge )
    {
        cerr << "Max p-value: " << pValueMax << endl;
    }
    
    if ( warningCount > 0 && ! parameters.reads )
    {
    	warnKmerSize(parameters, *this, lengthMax, lengthMaxName, randomChance, kMin, warningCount);
    }
    
    return 0;
}

void CommandTriangle::writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValueMax) const
{
	const Sketch & sketch = output->sketch;
	const Sketch::Reference & ref = sketch.getReference(output->index);
	
    if ( !edge )
    {
        cout << (comment ? ref.comment : ref.name);
    }
	
    for ( uint64_t i = 0; i < output->index; i++ )
    {
        const CommandDistance::CompareOutput::PairOutput * pair = &output->pairs[i];
        
        if ( edge )
        {
            const Sketch::Reference & qry = sketch.getReference(i);
            cout << (comment ? qry.comment : qry.name) << '\t'<< (comment ? ref.comment : ref.name) << '\t' << pair->distance << endl;// << '\t' << pair->pValue << '\t' << pair->numer << '/' << pair->denom << endl;
        }
	    else
	    {
	        cout << '\t' << pair->distance;
	    }
	    
	    if ( pair->pValue > pValueMax )
	    {
	        pValueMax = pair->pValue;
	    }
    }
    
    if ( !edge )
    {
        cout << endl;
    }
    
    delete output;
}

CommandTriangle::TriangleOutput * compare(CommandTriangle::TriangleInput * input)
{
    const Sketch & sketch = input->sketch;
    
    CommandTriangle::TriangleOutput * output = new CommandTriangle::TriangleOutput(input->sketch, input->index);
    
    uint64_t sketchSize = sketch.getMinHashesPerWindow();
    
    for ( uint64_t i = 0; i < input->index; i++ )
    {
    	compareSketches(&output->pairs[i], sketch.getReference(input->index), sketch.getReference(i), sketchSize, sketch.getKmerSize(), sketch.getKmerSpace(), -1., -1.);
    }
    
    return output;
}

} // namespace mash
