#include <iostream>

#include "CommandList.h"

using namespace::std;

CommandList::~CommandList()
{
	for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
	{
		delete i->second;
	}
}

void CommandList::addCommand(Command * command)
{
	commands[command->name] = command;
}

void CommandList::print()
{
	for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
	{
		cout << i->first << '\t' << i->second->description << endl;
	}
}

int CommandList::run(int argc, const char ** argv)
{
	return commands.at(argv[1])->run(argc - 2, argv + 2);
}
