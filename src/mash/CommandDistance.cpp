#include "CommandDistance.h"
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"
#include "gsl/gsl_cdf.h"

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
    useOption("bloomError");
    addOption("table", Option(Option::Boolean, "t", "Table output.", ""));
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
    bool table = options.at("table").active;
    
    Sketch::Parameters parameters;
    
    parameters.kmerSize = options.at("kmer").getArgumentAsNumber();
    parameters.minHashesPerWindow = options.at("sketchSize").getArgumentAsNumber();
    parameters.concatenated = options.at("concat").active;
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
        
        threadPool.runWhenThreadAvailable(new CompareInput(sketch, sketchQuery, queryFiles[i], parameters));
        
        while ( threadPool.outputAvailable() )
        {
            writeOutput(threadPool.popOutputWhenAvailable(), table);
        }
    }
    
    while ( threadPool.running() )
    {
        writeOutput(threadPool.popOutputWhenAvailable(), table);
    }
    
    return 0;
}

void CommandDistance::writeOutput(CompareOutput * output, bool table) const
{
    for ( int i = 0; i < output->pairs.size(); i++ )
    {
        const CompareOutput::PairOutput & pair = output->pairs.at(i);
        string queryLast;
        
        if ( table )
        {
            if ( i > 0 && pair.nameQuery != output->pairs.at(i - 1).nameQuery )
            {
                cout << endl << pair.nameQuery << '\t' << pair.score;
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
                
                cout << pair.score;
            }
        }
        else
        {
            cout << pair.score << '\t' << pair.nameRef << '\t' << pair.nameQuery << '\t' << pair.pValue << endl;
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
            
            compareSketches(output->pairs[pairIndex], sketchRef.getReference(j), sketchQuery->getReference(i), sketchSize, sketchRef.getKmerSpace());
        }
    }
    
    delete data->sketchQuery;
    
    return output;
}

void compareSketches(CommandDistance::CompareOutput::PairOutput & output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, int sketchSize, uint64_t kmerSpace)
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
    
    output.score = 1. - float(common) / denom;
    output.pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);
    
    output.nameRef = refRef.name;
    output.nameQuery = refQry.name;
}

double pValue(int x, int lengthRef, int lengthQuery, uint64_t kmerSpace, int sketchSize)
{
    double pX = 1. / (1. + (double)kmerSpace / lengthRef);
    double pY = 1. / (1. + (double)kmerSpace / lengthQuery);
    
    double r = pX * pY / (pX + pY - pX * pY);
    
    double M = (double)kmerSpace * (pX + pY) / (1. + r);
    
    return gsl_cdf_hypergeometric_Q(x - 1, r * M, M - r * M, sketchSize);
}
