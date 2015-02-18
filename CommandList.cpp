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
	cout << endl << "minimap";
	
	for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
	{
		cout << '\t' << i->first << '\t' << i->second->description << endl;
	}
	
	cout << endl;
}

int CommandList::run(int argc, const char ** argv)
{
	if ( argc < 2 || commands.count(argv[1]) == 0 )
	{
		print();
		return 0;
	}
	
	return commands.at(argv[1])->run(argc - 2, argv + 2);
}
