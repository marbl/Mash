// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

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
