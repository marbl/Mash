// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

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
