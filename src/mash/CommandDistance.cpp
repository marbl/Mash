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
#include <math.h>

#ifdef USE_BOOST
    #include <boost/math/distributions/binomial.hpp>
    using namespace::boost::math;
#else
    #include <gsl/gsl_cdf.h>
#endif

using namespace::std;

CommandDistance::CommandDistance()
: Command()
{
    name = "dist";
    summary = "Estimate the distance of query sequences to references.";
    description = "Estimate the distance of each query sequence to the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or Mash sketch files (.msh) with matching k-mer sizes. Query files can also be files of file names (see -l). Whole files are compared by default (see -i). The output fields are [reference-ID, query-ID, distance, p-value, shared-hashes].";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    useOption("threads");
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Each query file contains a list of sequence files, one per line. The reference file is not affected.", ""));
    addOption("table", Option(Option::Boolean, "t", "Output", "Table output (will not report p-values, but fields will be blank if they do not meet the p-value threshold).", ""));
    //addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances and divide by k-mer size to provide a better analog to phylogenetic distance. The special case of zero shared min-hashes will result in a distance of 1.", ""));
    addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
    addOption("distance", Option(Option::Number, "d", "Output", "Maximum distance to report.", "1.0", 0., 1.));
    useOption("kmer");
    useOption("sketchSize");
    useOption("individual");
    useOption("warning");
    useOption("noncanonical");
    useOption("reads");
    useOption("minCov");
}

int CommandDistance::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        /*char tab = '\t';
        
        for ( int kmerSize = 4; kmerSize <= 32; kmerSize++ )
        {
            double kmerSpace = pow(4, kmerSize);
            
            for ( uint64_t refSize = 10000; refSize <= 1000000000000; refSize *= 10 )
            {
                for ( uint64_t qrySize = 10000; qrySize <= refSize; qrySize *= 10 )
                {
                    for ( int sketchSize = 100; sketchSize <= 1000; sketchSize += 100 )
                    {
                        for ( int common = 1; common <= sketchSize + 1; common += 10 )
                        {
                            if ( common > sketchSize )
                            {
                                common = sketchSize;
                            }
                            
                            if ( common > kmerSpace )
                            {
                                continue;
                            }
                            
                            double pX = 1. / (1. + (double)kmerSpace / refSize);
                            double pY = 1. / (1. + (double)kmerSpace / qrySize);
    
                            double r = pX * pY / (pX + pY - pX * pY);
    
                            uint64_t M = (double)kmerSpace * (pX + pY) / (1. + r);
                            
                            //cout << "k: " << kmerSize << tab << "L1: " << refSize << tab << "L2: " << qrySize << tab << "s: " << sketchSize << tab << "x: " << common << tab << " | " << "Ek: " << kmerSpace << tab << "pX: " << pX << tab << "pY: " << pY << tab << "r: " << r << tab << "M: " << M << tab;
                            //cout << (M < sketchSize ? M : sketchSize) << tab << r * M << tab << M - r * M << endl;
                            //double p = cdf(complement(hypergeometric_distribution(r * M, M < sketchSize ? M : sketchSize, M), common - 1 ));
                            //double p = cdf(complement(binomial(M < sketchSize ? M : sketchSize, r), common - 1 ));
                            //double p = gsl_cdf_hypergeometric_Q(common - 1, r * M, M - uint64_t(r * M), M < sketchSize ? M : sketchSize);
                            double p = gsl_cdf_binomial_Q(common - 1, r, M < sketchSize ? M : sketchSize);
                            
                            cout << p << endl;
                        }
                    }
                }
            }
        }*/
        
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
    
    parameters.kmerSize = options.at("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = options.at("sketchSize").getArgumentAsNumber();
    parameters.concatenated = ! options.at("individual").active;
    parameters.noncanonical = options.at("noncanonical").active;
    parameters.reads = options.at("reads").active;
    parameters.minCov = options.at("minCov").getArgumentAsNumber();
    parameters.warning = options.at("warning").getArgumentAsNumber();
    
    if ( options.at("minCov").active )
    {
        parameters.reads = true;
    }
    
    if ( parameters.reads && ! parameters.concatenated )
    {
        cerr << "ERROR: The option " << options.at("individual").identifier << " cannot be used with " << options.at("unique").identifier << "." << endl;
        return 1;
    }
    
    Sketch sketch;
    
    uint64_t lengthThreshold = (parameters.warning * pow(parameters.protein ? 20 : 4, parameters.kmerSize)) / (1. - parameters.warning);
    uint64_t lengthMax;
    double randomChance;
    int kMin;
    string lengthMaxName;
    int warningCount = 0;
    
    const string & fileReference = arguments[0];
    
    if ( hasSuffix(fileReference, suffixSketch) )
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
        
        sketch.initFromCapnp(fileReference.c_str());
        
        if ( options.at("sketchSize").active )
        {
            if ( parameters.reads && parameters.minHashesPerWindow != sketch.getMinHashesPerWindow() )
            {
                cerr << "ERROR: The sketch size must match the reference when using a bloom filter (leave this option out to inherit from the reference sketch)." << endl;
                return 1;
            }
        }
        else
        {
            parameters.minHashesPerWindow = sketch.getMinHashesPerWindow();
        }
        
        parameters.kmerSize = sketch.getKmerSize();
        parameters.noncanonical = sketch.getNoncanonical();
    }
    else
    {
        bool sketchFileExists = false;//sketch.initHeaderFromBaseIfValid(fileReference, false);
        /*
        if
        (
            (options.at("kmer").active && parameters.kmerSize != sketch.getKmerSize())
        )
        {
            sketchFileExists = false;
        }
        */
        if ( sketchFileExists )
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
            cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...";
            
            sketch.initFromSequence(refArgVector, parameters);
            
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
            
            cerr << "done.\n";
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
    
    if ( table )
    {
        cout << "#query";
        
        for ( int i = 0; i < sketch.getReferenceCount(); i++ )
        {
            cout << '\t' << sketch.getReference(i).name;
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
    
    for ( int i = 0; i < queryFiles.size(); i++ )
    {
        // If the input is a sketch file, load in the main thread; otherwise,
        // leave it to the child. Either way, the child will delete.
        //
        Sketch * sketchQuery = new Sketch();
        bool isSketch = hasSuffix(queryFiles[i], suffixSketch);
        
        if ( isSketch )
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
        
        threadPool.runWhenThreadAvailable(new CompareInput(sketch, sketchQuery, queryFiles[i], parameters, distanceMax, pValueMax));
        /*
        if ( ! isSketch )
        {
            for ( int j = 0; j < sketchQuery->getReferenceCount(); j++ )
            {
                int length = sketchQuery->getReference(j).length;
                
                if ( length > lengthThreshold )
                {
                    if ( warningCount == 0 || length > lengthMax )
                    {
                        lengthMax = length;
                        lengthMaxName = sketchQuery->getReference(j).name;
                        randomChance = sketchQuery->getRandomKmerChance(j);
                        kMin = sketchQuery->getMinKmerSize(j);
                    }
                    
                    warningCount++;
                }
            }
        }
        */
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
    	sketch.warnKmerSize(lengthMax, lengthMaxName, randomChance, kMin, warningCount);
    }
    
    return 0;
}

void CommandDistance::writeOutput(CompareOutput * output, bool table) const
{
    uint64_t refCount = output->sketchRef.getReferenceCount();
    
    for ( uint64_t i = 0; i < output->sketchQuery->getReferenceCount(); i++ )
    {
        if ( table )
        {
            cout << output->sketchQuery->getReference(i).name;
        }
        
        for ( uint64_t j = 0; j < refCount; j++ )
        {
            const CompareOutput::PairOutput & pair = output->pairs.at(i * refCount + j);
        
            if ( table )
            {
                cout << '\t';
            
                if ( pair.pass )
                {
                    cout << pair.distance;
                }
            }
            else if ( pair.pass )
            {
                cout << output->sketchRef.getReference(j).name << '\t' << output->sketchQuery->getReference(i).name << '\t' << pair.distance << '\t' << pair.pValue << '\t' << pair.numer << '/' << pair.denom << endl;
            }
        }
    
        if ( table )
        {
            cout << endl;
        }
    }
    
    delete output;
}

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * data)
{
    const Sketch & sketchRef = data->sketchRef;
    Sketch * sketchQuery = data->sketchQuery;
    
    CommandDistance::CompareOutput * output = new CommandDistance::CompareOutput(data->sketchRef, data->sketchQuery);
    
    if ( sketchQuery->getReferenceCount() == 0 )
    {
        // input was sequence file; sketch now
        
        vector<string> fileVector;
        fileVector.push_back(data->file);
        
        sketchQuery->initFromSequence(fileVector, data->parameters);
    }
    
    uint64_t sketchSize = sketchQuery->getMinHashesPerWindow() < sketchRef.getMinHashesPerWindow() ?
        sketchQuery->getMinHashesPerWindow() :
        sketchRef.getMinHashesPerWindow();
    
    output->pairs.resize(sketchRef.getReferenceCount() * sketchQuery->getReferenceCount());
    
    for ( uint64_t i = 0; i < sketchQuery->getReferenceCount(); i++ )
    {
        for ( uint64_t j = 0; j < sketchRef.getReferenceCount(); j++ )
        {
            uint64_t pairIndex = i * sketchRef.getReferenceCount() + j;
            
            compareSketches(output->pairs[pairIndex], sketchRef.getReference(j), sketchQuery->getReference(i), sketchSize, sketchRef.getKmerSize(), sketchRef.getKmerSpace(), data->maxDistance, data->maxPValue);
        }
    }
    
    return output;
}

void compareSketches(CommandDistance::CompareOutput::PairOutput & output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxDistance, double maxPValue)
{
    uint64_t i = 0;
    uint64_t j = 0;
    uint64_t common = 0;
    uint64_t denom = 0;
    const HashList & hashesSortedRef = refRef.hashesSorted;
    const HashList & hashesSortedQry = refQry.hashesSorted;
    
    output.pass = false;
    
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
    
    if ( distance > maxDistance )
    {
        return;
    }
    
    output.numer = common;
    output.denom = denom;
    output.distance = distance;
    output.pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);
    
    if ( output.pValue > maxPValue )
    {
        return;
    }
    
    output.pass = true;
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
