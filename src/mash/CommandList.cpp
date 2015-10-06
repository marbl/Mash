#include <iostream>

#include "CommandList.h"
#include "version.h"

using namespace::std;

CommandList::CommandList(string nameNew)
{
    name = nameNew;
}

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
    vector<vector<string>> columns(1);
    
    cout << endl << "Version: " << version << endl;
    cout << endl << "Usage:" << endl << endl;
    
    columns[0].push_back(name + " <command> [options] [arguments ...]");
    printColumns(columns);
    
    cout << "Commands:" << endl << endl;
    
    int lengthMax = 0;
    
    for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
    {
        if ( i->first.length() > lengthMax )
        {
            lengthMax = i->first.length();
        }
    }
    
    columns.clear();
    columns.resize(2);
    
    for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
    {
        columns[0].push_back(i->first);
        columns[1].push_back(i->second->summary);
    }
    
    printColumns(columns);
}

int CommandList::run(int argc, const char ** argv)
{
	if ( argc > 1 && strcmp(argv[1], "--version") == 0 )
	{
		cout << version << endl;
		return 0;
	}
	
    if ( argc < 2 || commands.count(argv[1]) == 0 )
    {
        print();
        return 0;
    }
    
    return commands.at(argv[1])->run(argc - 2, argv + 2);
}
