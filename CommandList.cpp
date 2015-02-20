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
	cout << endl << "Usage:" << endl << endl;
	cout << "   minimap <command> [options] [arguments ...]" << endl << endl;
	cout << "Commands:" << endl << endl;
	
	int lengthMax = 0;
	
	for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
	{
		if ( i->first.length() > lengthMax )
		{
			lengthMax = i->first.length();
		}
	}
	
	vector<vector<string>> columns(2);
	
	for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
	{
		columns[0].push_back(i->first);
		columns[1].push_back(i->second->description);
	}
	
	printColumns(columns);
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
