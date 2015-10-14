// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

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
