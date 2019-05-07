// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandTriangle
#define INCLUDED_CommandTriangle

#include "Command.h"
#include "CommandDistance.h"
#include "Sketch.h"

namespace mash {

class CommandTriangle : public Command
{
public:
    
    struct TriangleInput
    {
        TriangleInput(const Sketch & sketchNew, uint64_t indexNew, const Sketch::Parameters & parametersNew)
            :
            sketch(sketchNew),
            index(indexNew),
            parameters(parametersNew)
            {}
        
        const Sketch & sketch;
        uint64_t index;
        const Sketch::Parameters & parameters;
    };
    
    struct TriangleOutput
    {
        TriangleOutput(const Sketch & sketchNew, uint64_t indexNew)
            :
            sketch(sketchNew),
            index(indexNew)
        {
            pairs = new CommandDistance::CompareOutput::PairOutput[index];
        }
        
        ~TriangleOutput()
        {
            delete [] pairs;
        }
        
        const Sketch & sketch;
        uint64_t index;
        
        CommandDistance::CompareOutput::PairOutput * pairs;
    };
    
    CommandTriangle();
    
    int run() const; // override
    
private:
    
    double pValueMax;
    bool comment;
    
    void writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValueMax) const;
};

CommandTriangle::TriangleOutput * compare(CommandTriangle::TriangleInput * input);

} // namespace mash

#endif
