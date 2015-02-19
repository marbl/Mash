#include "CommandList.h"
#include "CommandIndex.h"

using namespace::std;

int main(int argc, const char ** argv)
{
	CommandList commandList;
	
	commandList.addCommand(new CommandIndex());
	
	return commandList.run(argc, argv);
}
