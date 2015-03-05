#ifndef INCLUDED_CommandCompare
#define INCLUDED_CommandCompare

#include "Command.h"
#include "Index.h"

class CommandCompare : public Command
{
public:
    
    struct CompareInput
    {
        CompareInput(const Index & indexRefNew, const std::string fileNew)
            :
            indexRef(indexRefNew),
            file(fileNew)
            {}
        
        const Index & indexRef;
        const std::string file;
    };
    
    struct CompareOutput
    {
        float score;
        std::string file;
    };
    
    CommandCompare();
    
    int run() const; // override
    
private:
    
    void writeOutput(CompareOutput * output) const;
};

CommandCompare::CompareOutput * compare(CommandCompare::CompareInput * data);

#endif
