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
    //commandList.addCommand(new CommandContain());
    commandList.addCommand(new CommandInfo());
    commandList.addCommand(new CommandPaste());
    
    return commandList.run(argc, argv);
}
