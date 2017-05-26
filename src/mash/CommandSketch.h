// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandSketch
#define INCLUDED_CommandSketch

#include "Command.h"

namespace mash {

class CommandSketch : public Command
{
public:

    CommandSketch();
    
    int run() const; // override
};

} // namespace mash

#endif
