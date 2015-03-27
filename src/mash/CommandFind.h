#ifndef INCLUDED_CommandFind
#define INCLUDED_CommandFind

#include "Command.h"
#include "Index.h"
#include <string.h>
#include <queue>

class CommandFind : public Command
{
public:
    
    struct FindInput
    {
        FindInput
        (
            const Index & indexNew,
            const char * seqIdNew,
            const char * seqNew,
            uint32_t lengthNew,
            float thresholdNew,
            int bestNew,
            bool selfMatchesNew
        ) :
        index(indexNew),
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
        
        const Index & index;
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
    
    void writeOutput(const Index & index, FindOutput * output) const;
};

CommandFind::FindOutput * find(CommandFind::FindInput * data);
void findPerStrand(const CommandFind::FindInput * input, CommandFind::FindOutput * output, bool minusStrand);
bool operator<(const CommandFind::FindOutput::Hit & a, const CommandFind::FindOutput::Hit & b);

#endif
