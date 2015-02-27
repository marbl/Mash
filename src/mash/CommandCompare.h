#ifndef INCLUDED_CommandCompare
#define INCLUDED_CommandCompare

#include "Command.h"
#include "Index.h"

class CommandCompare : public Command
{
public:

    CommandCompare();
    
    int run() const; // override
};

float compare(const Index & indexRef, const std::string file);

#endif
