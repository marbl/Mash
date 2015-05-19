#include "CommandList.h"
#include "CommandSketch.h"
#include "CommandFind.h"
#include "CommandDistance.h"

using namespace::std;

int main(int argc, const char ** argv)
{
    CommandList commandList("mash");
    
    commandList.addCommand(new CommandSketch());
    //commandList.addCommand(new CommandFind());
    commandList.addCommand(new CommandDistance());
    
    return commandList.run(argc, argv);
}
