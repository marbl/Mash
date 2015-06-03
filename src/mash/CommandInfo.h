#ifndef INCLUDED_CommandInfo
#define INCLUDED_CommandInfo

#include "Command.h"
#include "Sketch.h"

class CommandInfo : public Command
{
public:
    
    CommandInfo();
    
    int run() const; // override
};

#endif
