#include "CommandList.h"
#include "CommandIndex.h"

using namespace::std;

int main(int argc, const char ** argv)
{
	CommandList commandList;
	
	Command * commandIndex = new CommandIndex("index", "Create a reference index");
	
	commandIndex->addOption("kmer", Command::Option(Command::Option::Number, "k", "Kmer size", "11"));
	commandIndex->addOption("factor", Command::Option(Command::Option::Number, "c", "Compression factor", "100"));
	commandIndex->addOption("hashes", Command::Option(Command::Option::Number, "m", "Number of hashes", "10"));
	
	commandList.addCommand(commandIndex);
	
	return commandList.run(argc, argv);
}
