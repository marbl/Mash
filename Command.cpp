#include <iostream>

#include "Command.h"

using namespace::std;

void Command::addOption(string name, Option option)
{
	options[name] = option;
	optionNamesByIdentifier[option.identifier] = name;
}

void Command::print()
{
	cout << endl << name;
	
	for ( map<string, Option>::iterator i = options.begin(); i != options.end(); i++ )
	{
		cout << "   -" << i->first << '\t';
		
		switch ( i->second.type )
		{
			case Option::Boolean:
				break;
			case Option::Number:
				cout << "<number>";
				break;
			case Option::File:
				cout << "<file>";
				break;
		}
		
		cout << '\t' << description << endl;
	}
}

int Command::run(int argc, const char ** argv)
{
	for ( int i = 0; i < argc; i++ )
	{
		if ( argv[i][0] == '-' )
		{
			Option & option = options.at(optionNamesByIdentifier.at(argv[i] + 1));
			
			if ( option.type != Option::Boolean )
			{
				i++;
				option.argument = argv[i];
			}
		}
		else
		{
			arguments.push_back(argv[i]);
		}
	}
	
	return run();
}
