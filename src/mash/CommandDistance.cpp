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
    description = "Estimate the distance of each query sequence (or file with -f) to the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or mash sketch files (.msh) with matching kmer sizes (-k). The distance is one minus the Jaccard score for the set of min-hashes whose size is that of the smaller sketch. The output fields are [distance, p-value, reference-ID, query-ID].";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    useOption("threads");
    addOption("list", Option(Option::Boolean, "l", "Input", "Query files are lists of file names.", ""));
    addOption("table", Option(Option::Boolean, "t", "Output", "Table output (will not report p-values, but fields will be blank if they do not meet the p-value threshold).", ""));
    addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances. To avoid taking the log of 0, 1-offset pseudocounts are used for the Jaccard index.", ""));
    addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
    useOption("kmer");
    useOption("sketchSize");
    useOption("individual");
    useOption("noncanonical");
    useOption("unique");
    useOption("genome");
    useOption("memory");
    useOption("bloomError");
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
    bool log = options.at("log").active;
    double pValueMax = options.at("pvalue").getArgumentAsNumber();
    
    Sketch::Parameters parameters;
    
    parameters.kmerSize = options.at("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = options.at("sketchSize").getArgumentAsNumber();
    parameters.concatenated = ! options.at("individual").active;
    parameters.noncanonical = options.at("noncanonical").active;
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
        
        threadPool.runWhenThreadAvailable(new CompareInput(sketch, sketchQuery, queryFiles[i], parameters, log));
        
        while ( threadPool.outputAvailable() )
        {
            writeOutput(threadPool.popOutputWhenAvailable(), table, pValueMax);
        }
    }
    
    while ( threadPool.running() )
    {
        writeOutput(threadPool.popOutputWhenAvailable(), table, pValueMax);
    }
    
    return 0;
}

void CommandDistance::writeOutput(CompareOutput * output, bool table, double pValueMax) const
{
    for ( int i = 0; i < output->pairs.size(); i++ )
    {
        const CompareOutput::PairOutput & pair = output->pairs.at(i);
        string queryLast;
        
        double score = pair.score;
        
        if ( table )
        {
            if ( i > 0 && pair.nameQuery != output->pairs.at(i - 1).nameQuery )
            {
                cout << endl << pair.nameQuery << '\t';
                
                if ( pair.pValue <= pValueMax )
                {
                    cout << score;
                }
            }
            else
            {
                if ( i == 0 )
                {
                    cout << pair.nameQuery << '\t';
                }
                else
                {
                    cout << '\t';
                }
                
                if ( pair.pValue <= pValueMax )
                {
                    cout << score;
                }
            }
        }
        else if ( pair.pValue <= pValueMax )
        {
            cout << score << '\t' << pair.pValue << '\t' << pair.nameRef << '\t' << pair.nameQuery << endl;
        }
    }
    
    if ( table )
    {
        cout << endl;
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
    
    int sketchSize = sketchQuery->getMinHashesPerWindow() < sketchRef.getMinHashesPerWindow() ?
        sketchQuery->getMinHashesPerWindow() :
        sketchRef.getMinHashesPerWindow();
    
    output->pairs.resize(sketchRef.getReferenceCount() * sketchQuery->getReferenceCount());
    
    for ( int i = 0; i < sketchQuery->getReferenceCount(); i++ )
    {
        for ( int j = 0; j < sketchRef.getReferenceCount(); j++ )
        {
            int pairIndex = i * sketchRef.getReferenceCount() + j;
            
            compareSketches(output->pairs[pairIndex], sketchRef.getReference(j), sketchQuery->getReference(i), sketchSize, sketchRef.getKmerSpace(), data->log);
        }
    }
    
    delete data->sketchQuery;
    
    return output;
}

void compareSketches(CommandDistance::CompareOutput::PairOutput & output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, int sketchSize, double kmerSpace, bool log)
{
    int i = 0;
    int j = 0;
    int common = 0;
    int denom = 0;
    
    const HashList & hashesSortedRef = refRef.hashesSorted;
    const HashList & hashesSortedQry = refQry.hashesSorted;
    
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
    
    int offset = log ? 1 : 0;
    
    if ( log )
    {
        output.score = -log10(double(common + 1) / (denom + 1));
    }
    else
    {
        output.score = 1. - double(common) / denom;
    }
    
    output.pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);
    
    output.nameRef = refRef.name;
    output.nameQuery = refQry.name;
}

double pValue(uint32_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint32_t sketchSize)
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
