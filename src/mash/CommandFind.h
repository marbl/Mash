#ifndef INCLUDED_CommandFind
#define INCLUDED_CommandFind

#include "Command.h"
#include "Index.h"
#include <string.h>

class CommandFind : public Command
{
public:
    
    struct FindInput
    {
        FindInput(const Index & indexNew, const char * seqIdNew, const char * seqNew, uint32_t lengthNew, float thresholdNew)
            :
            index(indexNew),
            length(lengthNew),
            threshold(thresholdNew),
            seqId(seqIdNew)
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
    };
    
    struct FindOutput
    {
        struct Hit
        {
            unsigned int ref;
            unsigned int start;
            unsigned int end;
            float score;
        };
        
        std::string seqId;
        std::vector<Hit> hits;
    };
    
    CommandFind();
    
    int run() const; // override
	
private:
	
	void writeOutput(const Index & index, const FindOutput * output) const;
};

CommandFind::FindOutput * find(CommandFind::FindInput * data);

#endif
