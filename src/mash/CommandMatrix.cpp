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
std::vector<double> compute_matrix(const Sketch& sketch, int threads);
void print_matrix(const std::vector<std::string>& names, const std::vector<double>& matrix);

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
    int bootstrap = 0;
    
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
    
    auto genomes = std::vector<Genome>{};
    std::transform(filenames.begin(), filenames.end(), std::back_inserter(genomes),[](std::string file_name){
        return Genome::fromFile(file_name);
    });

    auto sketchAll = Sketch(genomes, parameters);
    auto count = sketchAll.getReferenceCount();
    
    auto names = std::vector<std::string>();
    names.reserve(count);

    for (size_t i = 0; i < count; i++) {
        names.push_back(extract_name(sketchAll.getReference(i).name));
    }

    auto matrix = compute_matrix(sketchAll, threads);
    print_matrix(names, matrix);

    if (bootstrap > 0) {
        // resketch
        auto sketch = Sketch(genomes, parameters);
        auto matrix = compute_matrix(sketch, threads);
        print_matrix(names, matrix);
    }

    return 0;
}


void print_matrix(const std::vector<std::string>& names, const std::vector<double>& matrix)
{
    auto count = names.size();
    auto mat = [&](const size_t i, const size_t j)->const double & {
        return matrix[i * count + j];
    };

    fprintf(stdout, "%zu\n", count);

    for ( uint64_t i = 0; i < count; i++ )
    {
        fprintf(stdout, "%-10s", names[i].c_str());

        for ( uint64_t j = 0; j < count; j++ )
        {
            fprintf(stdout, " %1.4e", mat(i,j));
        }
        fprintf(stdout, "\n");
    }
}

std::vector<double> compute_matrix(const Sketch& sketch, int threads)
{
    auto ret = std::vector<double>();
    auto count = sketch.getReferenceCount();
    std::cerr << "reference count: " << count << std::endl;
    ret.reserve(count * count);
    auto mat = [&](const size_t i, const size_t j)->double & {
        return ret[i * count + j];
    };

    #pragma omp parallel for num_threads(threads)
    for ( uint64_t i = 0; i < count; i++ )
    {
        mat(i,i) = 0.0;

        for ( uint64_t j = 0; j < i; j++ )
        {
            // run
            auto output = compare(sketch, j, i);
            mat(j,i) = mat(i,j) = output.distance;
        }
    }

    return ret;
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
