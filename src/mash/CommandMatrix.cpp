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
};

PairOutput compare(const Sketch & sketchNew, uint64_t indexRefNew, uint64_t indexQueryNew);

CommandMatrix::CommandMatrix()
: Command()
{
    name = "matrix";
    summary = "Compute the distance matrix.";
    description = "Estimate the distances of sequences and produce a matrix (PHYLIP style).";
    argumentString = "<file> ...";
    
    useOption("help");
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <file> specify paths to sequence files, one per line.", ""));
    // addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
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
    
    auto seqs = std::vector<Sequence>{};
    std::transform(filenames.begin(), filenames.end(), std::back_inserter(seqs),[](std::string file_name){
        return Sequence::fromFile(file_name);
    });

    auto sketchAll = Sketch(seqs, parameters);
    
    uint64_t count = sketchAll.getReferenceCount();

    auto matbase = new double[count*count];
    auto mat = new double*[count];

    for ( uint64_t i = 0; i < count; i++ )
    {
        mat[i] = &matbase[i * count];
    }

    #pragma omp parallel for num_threads(threads)
    for ( uint64_t i = 0; i < count; i++ )
    {
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
    
    uint64_t i = 0;
    uint64_t j = 0;
    uint64_t common = 0;
    uint64_t denom = 0;
    const HashList & hashesSortedRef = sketch.getReference(indexRef).hashesSorted;
    const HashList & hashesSortedQry = sketch.getReference(indexQuery).hashesSorted;
    
    /*
     * The following code should be replaced with std::set_intersection, but
     * some poor design choices of HashList prevent this.
     */

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
    
    /*
     * If I am not mistaken, this block could be just:
     *
     *  denom = hashesSortedRef.size() + hashesSortedQry.size() - common;
     *
     * But I am not sure in the cases of different sizes.
     */

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
        distance = -log(2 * jaccard / (1. + jaccard)) / sketch.getKmerSize();
    }

    PairOutput output;
    
    output.numer = common;
    output.denom = denom;
    output.distance = distance;
    
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
