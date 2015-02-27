#ifndef INCLUDED_CommandIndex
#define INCLUDED_CommandIndex

#include "Command.h"

class CommandIndex : public Command
{
public:

    CommandIndex();
    
    int run() const; // override
};

#endif
