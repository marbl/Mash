#ifndef INCLUDED_CommandPaste
#define INCLUDED_CommandPaste

#include "Command.h"
#include "Sketch.h"

class CommandPaste : public Command
{
public:
    
    CommandPaste();
    
    int run() const; // override
};

#endif
