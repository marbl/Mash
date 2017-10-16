// Copyright © 2015,2017 Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, Adam Phillippy, and Fabian Klötzl
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandMatrix.h"
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

CommandMatrix::CommandMatrix()
: Command()
{
    name = "matrix";
    summary = "Compute the distance matrix.";
    description = "Estimate the distances of sequences and produce a matrix (PHYLIP style).";
    argumentString = "<file> ...";
    
    useOption("help");
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <query> specify paths to sequence files, one per line. The reference file is not affected.", ""));
    addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
    useSketchOptions();
}

int CommandMatrix::run() const
{
    if ( arguments.size() < 1 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    bool list = options.at("list").active;
    //bool log = options.at("log").active;
    double pValueMax = options.at("pvalue").getArgumentAsNumber();
    
    Sketch::Parameters parameters;
    
    if ( sketchParameterSetup(parameters, *(Command *)this) )
    {
    	return 1;
    }
    
    uint64_t lengthMax;
    double randomChance;
    int kMin;
    string lengthMaxName;
    int warningCount = 0;
    
    const string & fileReference = arguments[0];
    
    bool isSketch = hasSuffix(fileReference, suffixSketch);
    
    vector<string> filenames;

    if ( list )
    {
        for ( string arg: arguments )
        {
            splitFile(arg, filenames);
        }
    }
    else
    {
        filenames = arguments;
    }
    
    
    Sketch sketchAll;
    sketchAll.initFromFiles(filenames, parameters, 0, true);
       
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

        if ( options.at("sketchSize").active )
        {
            if ( parameters.reads && parameters.minHashesPerWindow != sketchAll.getMinHashesPerWindow() )
            {
                cerr << "ERROR: The sketch size must match the reference when using a bloom filter (leave this option out to inherit from the reference sketch)." << endl;
                return 1;
            }
        }
        
        parameters.minHashesPerWindow = sketchAll.getMinHashesPerWindow();
        parameters.kmerSize = sketchAll.getKmerSize();
        parameters.noncanonical = sketchAll.getNoncanonical();
        parameters.preserveCase = sketchAll.getPreserveCase();
        parameters.seed = sketchAll.getHashSeed();
        
        string alphabet;
        sketchAll.getAlphabetAsString(alphabet);
        setAlphabetFromString(parameters, alphabet.c_str());
    }
    else
    {
        cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...";

        double lengthThreshold = (parameters.warning * sketchAll.getKmerSpace()) / (1. - parameters.warning);
        for ( uint64_t i = 0; i < sketchAll.getReferenceCount(); i++ )
        {
            uint64_t length = sketchAll.getReference(i).length;
        
            if ( length <= lengthThreshold )
            {
                continue;
            }

            if ( warningCount == 0 || length > lengthMax )
            {
                lengthMax = length;
                lengthMaxName = sketchAll.getReference(i).name;
                randomChance = sketchAll.getRandomKmerChance(i);
                kMin = sketchAll.getMinKmerSize(i);
            }
        
            warningCount++;
        }
    
        cerr << "done.\n";
    }
    
    ThreadPool<CompareInput, CompareOutput> threadPool(compare, threads);
    uint64_t count = sketchAll.getReferenceCount();
    uint64_t pairCount = count * count;
    uint64_t pairsPerThread = pairCount / parameters.parallelism;
    static uint64_t maxPairsPerThread = 0x1000;

    pairsPerThread = std::max(pairsPerThread, 1lu);
    pairsPerThread = std::min(pairsPerThread, maxPairsPerThread);
    
    uint64_t iFloor = pairsPerThread / count;
    uint64_t iMod = pairsPerThread % count;

    auto matbase = new double[count*count];
    auto mat = new double*[count];

    for ( uint64_t i = 0; i < count; i++ )
    {
        mat[i] = &matbase[i * count];
        mat[i][i] = 0.0;

        for ( uint64_t j = 0; j < i; j++ )
        {
            auto task = CompareInput(sketchAll, sketchAll, j, i, /*pairsPerThread*/1, parameters, pValueMax);
            // run
            auto output = compare(&task);
            mat[j][i] = mat[i][j] = output->pairs[0].distance;
        }
    }

    cout << count << "\n";
    for ( uint64_t i = 0; i < count; i++ )
    {
        cout << sketchAll.getReference(i).name;

        for ( uint64_t j = 0; j < count; j++ )
        {
            cout << " " << mat[i][j];
        }
        cout << "\n";
    }

    delete[] mat;
    delete[] matbase;
    
    return 0;
}

CommandMatrix::CompareOutput * compare(CommandMatrix::CompareInput * input)
{
    const Sketch & sketchRef = input->sketchRef;
    const Sketch & sketchQuery = input->sketchQuery;
    
    CommandMatrix::CompareOutput * output = new CommandMatrix::CompareOutput(input->sketchRef, input->sketchQuery, input->indexRef, input->indexQuery, input->pairCount);
    
    uint64_t sketchSize = sketchQuery.getMinHashesPerWindow() < sketchRef.getMinHashesPerWindow() ?
        sketchQuery.getMinHashesPerWindow() :
        sketchRef.getMinHashesPerWindow();
    
    uint64_t i = input->indexQuery;
    uint64_t j = input->indexRef;
    
    for ( uint64_t k = 0; k < input->pairCount && i < sketchQuery.getReferenceCount(); k++ )
    {
        compareSketches(&output->pairs[k], sketchRef.getReference(j), sketchQuery.getReference(i), sketchSize, sketchRef.getKmerSize(), sketchRef.getKmerSpace(), input->maxPValue);
        
        j++;
        
        if ( j == sketchRef.getReferenceCount() )
        {
            j = 0;
            i++;
        }
    }
    
    return output;
}

void compareSketches(CommandMatrix::CompareOutput::PairOutput * output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxPValue)
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
    }
    
    output->numer = common;
    output->denom = denom;
    output->distance = distance;
    output->pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);
    
    if ( output->pValue > maxPValue )
    {
        return;
    }
    
    output->pass = true;
}

} // namespace mash
