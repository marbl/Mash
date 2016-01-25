// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandContain
#define INCLUDED_CommandContain

#include "Command.h"
#include "Sketch.h"

class CommandContain : public Command
{
public:
    
    struct ContainInput
    {
        ContainInput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew, const Sketch::Parameters & parametersNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew),
            parameters(parametersNew)
            {}
        
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
		
        uint64_t indexRef;
        uint64_t indexQuery;
        
        std::string nameRef;
        const Sketch::Parameters & parameters;
    };
    
    struct ContainOutput
    {
        ContainOutput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew)
        {}
        
		float score;
		float error;
    
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
		
        uint64_t indexRef;
        uint64_t indexQuery;
    };
    
    CommandContain();
    
    int run() const; // override
    
private:
    
    void writeOutput(ContainOutput * output, float error) const;
};

CommandContain::ContainOutput * contain(CommandContain::ContainInput * data);
float containSketches(const HashList & hashesSortedRef, const HashList & hashesSortedQuery, float & errorToSet);

#endif
