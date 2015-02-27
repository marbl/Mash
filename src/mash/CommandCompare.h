#ifndef INCLUDED_CommandCompare
#define INCLUDED_CommandCompare

#include "Command.h"
#include "Index.h"

class CommandCompare : public Command
{
public:
    
    struct CompareData
    {
        CompareData(const Index & indexRefNew, const std::string fileNew)
            :
            indexRef(indexRefNew),
            file(fileNew)
            {}
        
        const Index & indexRef;
        const std::string file;
    };
    
    CommandCompare();
    
    int run() const; // override
};

void * compare(void *);

#endif
