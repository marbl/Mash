// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandList.h"
#include "CommandSketch.h"
#include "CommandFind.h"
#include "CommandDistance.h"
#include "CommandContain.h"
#include "CommandInfo.h"
#include "CommandPaste.h"

using namespace::std;

int main(int argc, const char ** argv)
{
    CommandList commandList("mash");
    
    commandList.addCommand(new CommandSketch());
    //commandList.addCommand(new CommandFind());
    commandList.addCommand(new CommandDistance());
#ifdef COMMAND_WITHIN
    commandList.addCommand(new CommandContain());
#endif
    commandList.addCommand(new CommandInfo());
    commandList.addCommand(new CommandPaste());
    
    return commandList.run(argc, argv);
}
