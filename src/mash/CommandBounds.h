// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandBounds
#define INCLUDED_CommandBounds

#include "Command.h"

namespace mash {

class CommandBounds : public Command
{
public:
    
    CommandBounds();
    int run() const; // override
};

} // namespace mash

#endif
