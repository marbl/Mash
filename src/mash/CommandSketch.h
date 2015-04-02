#ifndef INCLUDED_CommandSketch
#define INCLUDED_CommandSketch

#include "Command.h"

class CommandSketch : public Command
{
public:

    CommandSketch();
    
    int run() const; // override
};

#endif
