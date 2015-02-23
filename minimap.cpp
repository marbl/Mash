#include "CommandList.h"
#include "CommandIndex.h"
#include "CommandFind.h"

using namespace::std;

int main(int argc, const char ** argv)
{
    CommandList commandList;
    
    commandList.addCommand(new CommandIndex());
    commandList.addCommand(new CommandFind());
    
    return commandList.run(argc, argv);
}
