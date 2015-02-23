#ifndef INCLUDED_CommandFind
#define INCLUDED_CommandFind

#include "Command.h"

class CommandFind : public Command
{
public:

    CommandFind();
    
    int run() const; // override
};

#endif
