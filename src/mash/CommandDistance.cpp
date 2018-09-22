// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandDistance.h"
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

CommandDistance::CommandDistance()
: Command()
{
    name = "dist";
    summary = "Estimate the distance of query sequences to references.";
    description = "Estimate the distance of each query sequence to the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or Mash sketch files (.msh) with matching k-mer sizes. Query files can also be files of file names (see -l). Whole files are compared by default (see -i). The output fields are [reference-ID, query-ID, distance, p-value, shared-hashes].";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <query> specify paths to sequence files, one per line. The reference file is not affected.", ""));
    addOption("table", Option(Option::Boolean, "t", "Output", "Table output (will not report p-values, but fields will be blank if they do not meet the p-value threshold).", ""));
    //addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances and divide by k-mer size to provide a better analog to phylogenetic distance. The special case of zero shared min-hashes will result in a distance of 1.", ""));
    addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
    addOption("distance", Option(Option::Number, "d", "Output", "Maximum distance to report.", "1.0", 0., 1.));
    useSketchOptions();
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
    bool table = options.at("table").active;
    //bool log = options.at("log").active;
    double pValueMax = options.at("pvalue").getArgumentAsNumber();
    double distanceMax = options.at("distance").getArgumentAsNumber();
    
    Sketch::Parameters parameters;
    
    if ( sketchParameterSetup(parameters, *(Command *)this) )
    {
    	return 1;
    }
    
    Sketch sketchRef;
    
    uint64_t lengthMax;
    double randomChance;
    int kMin;
    string lengthMaxName;
    int warningCount = 0;
    
    const string & fileReference = arguments[0];
    
    bool isSketch = hasSuffix(fileReference, suffixSketch);
    
    if ( isSketch )
    {
        if ( options.at("kmer").active )
        {
            cerr << "ERROR: The option -" << options.at("kmer").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
            return 1;
        }
        
        if ( options.at("noncanonical").active )
        {
            cerr << "ERROR: The option -" << options.at("noncanonical").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
            return 1;
        }
        
        if ( options.at("protein").active )
        {
            cerr << "ERROR: The option -" << options.at("protein").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
            return 1;
        }
        
        if ( options.at("alphabet").active )
        {
            cerr << "ERROR: The option -" << options.at("alphabet").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
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
    
    double lengthThreshold = (parameters.warning * sketchRef.getKmerSpace()) / (1. - parameters.warning);
    
    if ( isSketch )
    {
        if ( options.at("sketchSize").active )
        {
            if ( parameters.reads && parameters.minHashesPerWindow != sketchRef.getMinHashesPerWindow() )
            {
                cerr << "ERROR: The sketch size must match the reference when using a bloom filter (leave this option out to inherit from the reference sketch)." << endl;
                return 1;
            }
        }
        
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
        for ( uint64_t i = 0; i < sketchRef.getReferenceCount(); i++ )
        {
            uint64_t length = sketchRef.getReference(i).length;
        
            if ( length > lengthThreshold )
            {
                if ( warningCount == 0 || length > lengthMax )
                {
                    lengthMax = length;
                    lengthMaxName = sketchRef.getReference(i).name;
                    randomChance = sketchRef.getRandomKmerChance(i);
                    kMin = sketchRef.getMinKmerSize(i);
                }
            
                warningCount++;
            }
        }
    
        cerr << "done.\n";
    }
    
    if ( table )
    {
        cout << "#query";
        
        for ( int i = 0; i < sketchRef.getReferenceCount(); i++ )
        {
            cout << '\t' << sketchRef.getReference(i).name;
        }
        
        cout << endl;
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
    
    Sketch sketchQuery;
    
    sketchQuery.initFromFiles(queryFiles, parameters, 0, true);
    
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
        
        threadPool.runWhenThreadAvailable(new CompareInput(sketchRef, sketchQuery, j, i, pairsPerThread, parameters, distanceMax, pValueMax));
        
        while ( threadPool.outputAvailable() )
        {
            writeOutput(threadPool.popOutputWhenAvailable(), table);
        }
    }
    
    while ( threadPool.running() )
    {
        writeOutput(threadPool.popOutputWhenAvailable(), table);
    }
    
    if ( warningCount > 0 && ! parameters.reads )
    {
    	warnKmerSize(parameters, *this, lengthMax, lengthMaxName, randomChance, kMin, warningCount);
    }
    
    return 0;
}

void CommandDistance::writeOutput(CompareOutput * output, bool table) const
{
    uint64_t i = output->indexQuery;
    uint64_t j = output->indexRef;
    
    for ( uint64_t k = 0; k < output->pairCount && i < output->sketchQuery.getReferenceCount(); k++ )
    {
        const CompareOutput::PairOutput * pair = &output->pairs[k];
        
        if ( table && j == 0 )
        {
            cout << output->sketchQuery.getReference(i).name;
        }
        
        if ( table )
        {
            cout << '\t';
    
            if ( pair->pass )
            {
                cout << pair->distance;
            }
        }
        else if ( pair->pass )
        {
            cout << output->sketchRef.getReference(j).name << '\t' << output->sketchQuery.getReference(i).name << '\t' << pair->distance << '\t' << pair->pValue << '\t' << pair->numer << '/' << pair->denom << endl;
        }
    
        j++;
        
        if ( j == output->sketchRef.getReferenceCount() )
        {
            if ( table )
            {
                cout << endl;
            }
            
            j = 0;
            i++;
        }
    }
    
    delete output;
}

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * input)
{
    const Sketch & sketchRef = input->sketchRef;
    const Sketch & sketchQuery = input->sketchQuery;
    
    CommandDistance::CompareOutput * output = new CommandDistance::CompareOutput(input->sketchRef, input->sketchQuery, input->indexRef, input->indexQuery, input->pairCount);
    
    uint64_t sketchSize = sketchQuery.getMinHashesPerWindow() < sketchRef.getMinHashesPerWindow() ?
        sketchQuery.getMinHashesPerWindow() :
        sketchRef.getMinHashesPerWindow();
    
    uint64_t i = input->indexQuery;
    uint64_t j = input->indexRef;
    
    for ( uint64_t k = 0; k < input->pairCount && i < sketchQuery.getReferenceCount(); k++ )
    {
        compareSketches(&output->pairs[k], sketchRef.getReference(j), sketchQuery.getReference(i), sketchSize, sketchRef.getKmerSize(), sketchRef.getKmerSpace(), input->maxDistance, input->maxPValue);
        
        j++;
        
        if ( j == sketchRef.getReferenceCount() )
        {
            j = 0;
            i++;
        }
    }
    
    return output;
}

void compareSketches(CommandDistance::CompareOutput::PairOutput * output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxDistance, double maxPValue)
{
    uint64_t i = 0;
    uint64_t j = 0;
    uint64_t common = 0;
    uint64_t denom = 0;
    const HashList & hashesSortedRef = refRef.hashesSorted;
    const HashList & hashesSortedQry = refQry.hashesSorted;
    
    output->pass = false;
    
    while ( denom < sketchSize && i < hashesSortedRef.size() && j < hashesSortedQry.size() )
    {
        if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQry.at(j), hashesSortedRef.get64()) )
        {
            i++;
        }
        else if ( hashLessThan(hashesSortedQry.at(j), hashesSortedRef.at(i), hashesSortedRef.get64()) )
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
    
    if ( denom < sketchSize )
    {
        // complete the union operation if possible
        
        if ( i < hashesSortedRef.size() )
        {
            denom += hashesSortedRef.size() - i;
        }
        
        if ( j < hashesSortedQry.size() )
        {
            denom += hashesSortedQry.size() - j;
        }
        
        if ( denom > sketchSize )
        {
            denom = sketchSize;
        }
    }
    
    double distance;
    double jaccard = double(common) / denom;
    
    if ( common == denom ) // avoid -0
    {
        distance = 0;
    }
    else if ( common == 0 ) // avoid inf
    {
        distance = 1.;
    }
    else
    {
        //distance = log(double(common + 1) / (denom + 1)) / log(1. / (denom + 1));
        distance = -log(2 * jaccard / (1. + jaccard)) / kmerSize;
        
        if ( distance > 1 )
        {
        	distance = 1;
        }
    }
    
    if ( maxDistance >= 0 && distance > maxDistance )
    {
        return;
    }
    
    output->numer = common;
    output->denom = denom;
    output->distance = distance;
    output->pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);
    
    if ( maxPValue >= 0 && output->pValue > maxPValue )
    {
        return;
    }
    
    output->pass = true;
}

double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize)
{
    if ( x == 0 )
    {
        return 1.;
    }
    
    double pX = 1. / (1. + kmerSpace / lengthRef);
    double pY = 1. / (1. + kmerSpace / lengthQuery);
    
    double r = pX * pY / (pX + pY - pX * pY);
    
    //double M = (double)kmerSpace * (pX + pY) / (1. + r);
    
    //return gsl_cdf_hypergeometric_Q(x - 1, r * M, M - r * M, sketchSize);
    
#ifdef USE_BOOST
    return cdf(complement(binomial(sketchSize, r), x - 1));
#else
    return gsl_cdf_binomial_Q(x - 1, r, sketchSize);
#endif
}

} // namespace mash
