#ifndef INCLUDED_CommandFind
#define INCLUDED_CommandFind

#include "Command.h"
#include "Index.h"

class CommandFind : public Command
{
public:

    CommandFind();
    
    int run() const; // override
};

void find(const Index & index, char * seq, uint32_t length, float threshold);

#endif
