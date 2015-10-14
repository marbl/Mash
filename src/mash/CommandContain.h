// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
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
        ContainInput(const Sketch & sketchRefNew, Sketch * sketchQueryNew, const std::string fileNew, const Sketch::Parameters & parametersNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            file(fileNew),
            parameters(parametersNew)
            {}
        
        const Sketch & sketchRef;
        Sketch * sketchQuery;
        std::string nameRef;
        const std::string file;
        const Sketch::Parameters & parameters;
    };
    
    struct ContainOutput
    {
        struct PairOutput
        {
            float score;
            float error;
            std::string nameRef;
            std::string nameQuery;
        };
        
        std::vector<PairOutput> pairs;
    };
    
    CommandContain();
    
    int run() const; // override
    
private:
    
    void writeOutput(ContainOutput * output, float error) const;
};

CommandContain::ContainOutput * contain(CommandContain::ContainInput * data);
float containSketches(const HashList & hashesSortedRef, const HashList & hashesSortedQuery, float & errorToSet);

#endif
