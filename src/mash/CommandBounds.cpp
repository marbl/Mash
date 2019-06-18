// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandBounds.h"
#include <iostream>
#include <math.h>

#ifdef USE_BOOST
    #include <boost/math/distributions/binomial.hpp>
    using namespace::boost::math;
#else
    #include <gsl/gsl_cdf.h>
#endif

using std::cout;
using std::endl;

namespace mash {

CommandBounds::CommandBounds()
: Command()
{
    name = "bounds";
    summary = "Print a table of Mash error bounds.";
    description = "Print a table of Mash error bounds for various sketch sizes and Mash distances based on a given k-mer size and desired confidence. Note that these calculations assume sequences are much larger than the sketch size, and will overestimate error bounds if this is not the case.";
    argumentString = "";
    
    useOption("help");
    addOption("kmer", Option(Option::Integer, "k", "", "k-mer size.", "21", 1, 32));
    addOption("prob", Option(Option::Number, "p", "", "Mash distance estimates will be within the given error bounds with this probability.", "0.99", 0, 1));
}

int CommandBounds::run() const
{
    if ( options.at("help").active )
    {
        print();
        return 0;
    }
    
	const int sketchSizeCount = 9;
	const double sketchSizes[] = {100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
	
	const int distCount = 8;
	const double dists[] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
	
	int k = getOption("kmer").getArgumentAsNumber();
	double q2 = (1.0 - getOption("prob").getArgumentAsNumber()) / 2.0;
	
	cout << endl << "Parameters (run with -h for details):" << endl;
	cout << "   k:   " << k << endl;
	cout << "   p:   " << getOption("prob").getArgumentAsNumber() << endl << endl;
	
	for ( int cont = 0; cont < 2; cont++ )
	{
		if ( cont )
		{
			cout << "\tScreen distance" << endl;
		}
		else
		{
			cout << "\tMash distance" << endl;
		}
		
		cout << "Sketch";
	
		for ( int i = 0; i < distCount; i++ )
		{
			cout << '\t' << dists[i];
		}
	
		cout << endl;
	
		for ( int i = 0; i < sketchSizeCount; i++ )
		{
			int s = sketchSizes[i];
			cout << s;
		
			for ( int j = 0; j < distCount; j++ )
			{
				double m2j;
				
				if ( cont )
				{
					//m2j = exp(-k * dists[j]);
					m2j = pow(1.0 - dists[j], k); // binomial model
				}
				else
				{
					m2j = 1.0 / (2.0 * exp(k * dists[j]) - 1.0);
				}
			
				int x = 0;
			
				while ( x < s )
				{
#ifdef USE_BOOST
					double cdfx = cdf(binomial(s, m2j), x);
#else
					double cdfx = gsl_cdf_binomial_P(x, m2j, s);
#endif
					if ( cdfx > q2 )
					{
						break;
					}
				
					x++;
				}
			
				double je = double(x) / s;
				double j2m;
				
				if ( cont )
				{
					//j2m = -1.0 / k * log(je);
					j2m = 1.0 - pow(je, 1. / k);
				}
				else
				{
					j2m = -1.0 / k * log(2.0 * je / (1.0 + je));
				}
				
				cout << '\t' << j2m - dists[j];
			}
		
			cout << endl;
		}
	
		cout << endl;
	}
	
	return 0;
}

} // namespace mash
