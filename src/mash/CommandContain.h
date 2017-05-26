// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandContain
#define INCLUDED_CommandContain

#include "Command.h"
#include "Sketch.h"

namespace mash {

class CommandContain : public Command
{
public:
    
    struct ContainInput
    {
        ContainInput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew, uint64_t pairCountNew, const Sketch::Parameters & parametersNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew),
            pairCount(pairCountNew),
            parameters(parametersNew)
            {}
        
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
		
        uint64_t indexRef;
        uint64_t indexQuery;
        uint64_t pairCount;
        
        std::string nameRef;
        const Sketch::Parameters & parameters;
    };
    
    struct ContainOutput
    {
        ContainOutput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew, uint64_t pairCountNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew),
            pairCount(pairCountNew)
        {
            pairs = new PairOutput[pairCount];
        }
        
        ~ContainOutput()
        {
            delete [] pairs;
        }
        
        struct PairOutput
        {
			double score;
			double error;
		};
    
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
		
        uint64_t indexRef;
        uint64_t indexQuery;
        uint64_t pairCount;
        
        PairOutput * pairs;
    };
    
    CommandContain();
    
    int run() const; // override
    
private:
    
    void writeOutput(ContainOutput * output, float error) const;
};

CommandContain::ContainOutput * contain(CommandContain::ContainInput * data);
double containSketches(const HashList & hashesSortedRef, const HashList & hashesSortedQuery, double & errorToSet);

} // namespace mash

#endif
