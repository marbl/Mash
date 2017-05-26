// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandFind
#define INCLUDED_CommandFind

#include "Command.h"
#include "Sketch.h"
#include <string.h>
#include <queue>

namespace mash {

class CommandFind : public Command
{
public:
    
    struct FindInput
    {
        FindInput
        (
            const Sketch & sketchNew,
            const char * seqIdNew,
            const char * seqNew,
            uint32_t lengthNew,
            float thresholdNew,
            int bestNew,
            bool selfMatchesNew
        ) :
        sketch(sketchNew),
        length(lengthNew),
        threshold(thresholdNew),
        seqId(seqIdNew),
        best(bestNew),
        selfMatches(selfMatchesNew)
        {
            seq = new char[strlen(seqNew) + 1];
            strcpy(seq, seqNew);
        }
        
        ~FindInput()
        {
            delete [] seq;
        }
        
        const Sketch & sketch;
        std::string seqId;
        char * seq;
        uint32_t length;
        float threshold;
        int best;
        bool selfMatches;
    };
    
    struct FindOutput
    {
        struct Hit
        {
            Hit(uint32_t refNew, uint32_t startNew, uint32_t endNew, bool minusStrandNew, float scoreNew)
                :
                ref(refNew),
                start(startNew),
                end(endNew),
                minusStrand(minusStrandNew),
                score(scoreNew)
                {}
            
            unsigned int ref;
            unsigned int start;
            unsigned int end;
            bool minusStrand;
            float score;
        };
        
        std::string seqId;
        std::priority_queue<Hit> hits;
    };
    
    CommandFind();
    
    int run() const; // override
    
private:
    
    void writeOutput(const Sketch & sketch, FindOutput * output) const;
};

CommandFind::FindOutput * find(CommandFind::FindInput * data);
void findPerStrand(const CommandFind::FindInput * input, CommandFind::FindOutput * output, bool minusStrand);
bool operator<(const CommandFind::FindOutput::Hit & a, const CommandFind::FindOutput::Hit & b);

} // namespace mash

#endif
