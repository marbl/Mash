// Copyright © 2015,2017 Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, Adam Phillippy, and Fabian Klötzl
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandMatrix.h"
#include "Sketch.h"
#include <iostream>
#include <vector>
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

std::string extract_name(const std::string &file_name);

struct PairOutput
{
    uint64_t numer;
    uint64_t denom;
    double distance;
    double pValue;
};

PairOutput compare(const Sketch & sketchNew, uint64_t indexRefNew, uint64_t indexQueryNew);
double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize);

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
    
    Sketch::Parameters parameters;
    
    if ( sketchParameterSetup(parameters, *(Command *)this) )
    {
        return 1;
    }
    
    vector<string> filenames;

    if ( options.at("list").active )
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
       
    bool isSketch = hasSuffix(filenames[0], suffixSketch);
    if ( isSketch )
    {
        auto unavailable_options = {"kmer", "noncanonical", "protein", "alphabet"};
        for ( auto opt: unavailable_options)
        {
            if ( options.at(opt).active )
            {
                cerr << "ERROR: The option -" << options.at(opt).identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
                return 1;
            }
        }

        if ( options.at("sketchSize").active )
        {
            if ( parameters.reads && parameters.minHashesPerWindow != sketchAll.getMinHashesPerWindow() )
            {
                cerr << "ERROR: The sketch size must match the reference when using a bloom filter (leave this option out to inherit from the reference sketch)." << endl;
                return 1;
            }
        }

        // As far as I can tell, `parameters` is unused.
    }
    
    uint64_t count = sketchAll.getReferenceCount();

    auto matbase = new double[count*count];
    auto mat = new double*[count];

    for ( uint64_t i = 0; i < count; i++ )
    {
        mat[i] = &matbase[i * count];
        mat[i][i] = 0.0;

        for ( uint64_t j = 0; j < i; j++ )
        {
            // run
            auto output = compare(sketchAll, j, i);
            mat[j][i] = mat[i][j] = output.distance;
        }
    }

    cout << count << "\n";
    for ( uint64_t i = 0; i < count; i++ )
    {
        auto file_name = sketchAll.getReference(i).name;
        cout << extract_name(file_name);

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

PairOutput compare(const Sketch & sketch, uint64_t indexRef, uint64_t indexQuery)
{
    uint64_t sketchSize = sketch.getMinHashesPerWindow();
    
    auto refRef = sketch.getReference(indexRef);
    auto refQry = sketch.getReference(indexQuery);
    int kmerSize = sketch.getKmerSize();
    double kmerSpace = sketch.getKmerSpace();

    uint64_t i = 0;
    uint64_t j = 0;
    uint64_t common = 0;
    uint64_t denom = 0;
    const HashList & hashesSortedRef = refRef.hashesSorted;
    const HashList & hashesSortedQry = refQry.hashesSorted;
    
    auto use64 = hashesSortedRef.get64();
    while ( denom < sketchSize && i < hashesSortedRef.size() && j < hashesSortedQry.size() )
    {
        if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQry.at(j), use64) )
        {
            i++;
        }
        else if ( hashLessThan(hashesSortedQry.at(j), hashesSortedRef.at(i), use64) )
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

    PairOutput output;
    
    output.numer = common;
    output.denom = denom;
    output.distance = distance;
    output.pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);
    
    return output;
}

std::string extract_name(const std::string &file_name)
{
    // find the last path separator
    auto left = file_name.rfind('/');
    left = (left == std::string::npos) ? 0 : left + 1;
    // left is the position one of to the right of the path separator

    // find the extension
    auto right = file_name.find('.', left);
    right = (right == std::string::npos) ? file_name.size() : right;

    // copy only the file name, not its path or extension
    return file_name.substr(left, right - left);
}

} // namespace mash
