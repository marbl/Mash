#include "CommandList.h"
#include "CommandIndex.h"
#include "CommandFind.h"
#include "CommandCompare.h"

using namespace::std;

int main(int argc, const char ** argv)
{
    CommandList commandList;
    
    commandList.addCommand(new CommandIndex());
    commandList.addCommand(new CommandFind());
    commandList.addCommand(new CommandCompare());
    
    return commandList.run(argc, argv);
}
