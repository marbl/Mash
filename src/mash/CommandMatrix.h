// Copyright © 2015,2017 Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, Adam Phillippy, and Fabian Klötzl
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandMatrix
#define INCLUDED_CommandMatrix

#include "Command.h"
#include "Sketch.h"

namespace mash {

class CommandMatrix : public Command
{
public:
    CommandMatrix();
    
    int run() const; // override    
};

} // namespace mash

#endif
