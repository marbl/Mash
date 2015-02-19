#include <iostream>

#include "Command.h"

using namespace::std;

float Command::Option::getArgumentAsNumber(float min, float max) const
{
	float number;
	bool failed = false;
	
	try
	{
		number = stof(argument);
		
		if ( min != max && (number < min || number > max) )
		{
			failed = true;
		}
	}
	catch ( const exception & e )
	{
		failed = true;
	}
	
	if ( failed )
	{
		cerr << "ERROR: Argument to -" << identifier << " must be a number";
		
		if ( min != max )
		{
			cerr << " between " << min << " and " << max;
		}
		
		cerr << " (" << argument << " given)" << endl;
		exit(1);
	}
	
	return number;
}

void Command::addOption(string name, Option option)
{
	options[name] = option;
	optionNamesByIdentifier[option.identifier] = name;
}

void Command::print() const
{
	cout << endl << "minimap " << name << " [options] " << argumentString << endl << endl;
	cout << description << endl << endl;
	cout << "Options:" << endl << endl;
	
	for ( map<string, Option>::const_iterator i = options.begin(); i != options.end(); i++ )
	{
		cout << "-" << i->second.identifier << '\t';
		
		switch ( i->second.type )
		{
			case Option::Boolean:
				break;
			case Option::Number:
				cout << "<number>";
				break;
			case Option::File:
				cout << "<path>  ";
				break;
		}
		
		cout << '\t' << i->second.description;
		
		if ( i->second.type != Option::Boolean )
		{
			cout << " [" << i->second.argument << ']';
		}
		
		cout << endl;
	}
	
	cout << endl;
}

int Command::run(int argc, const char ** argv)
{
	for ( int i = 0; i < argc; i++ )
	{
		if ( argv[i][0] == '-' )
		{
			if ( optionNamesByIdentifier.count(argv[i] + 1) == 0 )
			{
				cerr << "ERROR: Unrecognized option: " << argv[i] << endl;
				return 1;
			}
			
			Option & option = options.at(optionNamesByIdentifier.at(argv[i] + 1));
			
			if ( option.type != Option::Boolean )
			{
				i++;
				
				if ( i == argc )
				{
					cerr << "ERROR: -" << option.identifier << " requires an argument" << endl;
					return 1;
				}
				
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
