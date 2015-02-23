#include <iostream>
#include <sys/ioctl.h>

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
    cout << endl << "Usage:" << endl << endl;
    cout << "   minimap " << name << " [options] " << argumentString << endl << endl;
    cout << "Description:" << endl << endl;
    cout << "   " << description << endl << endl;
    
    if ( options.size() == 0 )
    {
        return;
    }
    
    cout << "Options:" << endl << endl;
    
    vector<vector<string>> columns(4);
    
    columns[0].push_back("Option");
    columns[1].push_back("Argument");
    columns[2].push_back("Default");
    columns[3].push_back("Description");
    
    for ( map<string, Option>::const_iterator i = options.begin(); i != options.end(); i++ )
    {
        columns[0].push_back("-" + i->second.identifier);
        
        string type;
        
        switch ( i->second.type )
        {
            case Option::Boolean:
                break;
            case Option::Number:
                type = "<number>";
                break;
            case Option::File:
                type = "<path>  ";
                break;
        }
        
        columns[1].push_back(type);
        
        string defaultString;
        
        if ( i->second.type != Option::Boolean )
        {
            defaultString = "[" + i->second.argument + "]";
        }
        
        columns[2].push_back(defaultString);
        columns[3].push_back(i->second.description);
    }
    
    printColumns(columns);
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
            
            if ( option.type == Option::Boolean )
            {
                option.active = true;
            }
            else
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

void printColumns(vector<vector<string>> columns, int indent, int spacing)
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    
    vector<int> lengthMaxes(columns.size(), 0);
    
    for ( int i = 0; i < columns.size(); i++ )
    {
        for ( int j = 0; j < columns[i].size(); j++ )
        {
            if ( columns[i][j].length() > lengthMaxes[i] )
            {
                lengthMaxes[i] = columns[i][j].length();
            }
        }
    }
    
    for ( int i = 0; i < columns[0].size(); i++ )
    {
        int offset = 0;
        int offsetTarget = indent;
        
        for ( int j = 0; j < columns.size(); j++ )
        {
            for ( int k = offset; k < offsetTarget; k++ )
            {
                cout << ' ';
            }
            
            int index = 0;
            const string & text = columns[j][i];
            
            do
            {
                int length = text.length() - index;
                
                if ( length + offsetTarget > w.ws_col )
                {
                    length = w.ws_col - offsetTarget;
                
                    while ( text[index + length] != ' ' && length > 0 )
                    {
                        length--;
                    }
                }
                
                if ( length == 0 )
                {
                    length = w.ws_col - offsetTarget;
                }
                
                if ( index > 0 )
                {
                    cout << endl;
                    
                    for ( int k = 0; k < offsetTarget; k++ )
                    {
                        cout << ' ';
                    }
                }
                
                cout << text.substr(index, length);
                index += length;
                
                while ( index < text.length() && text[index] == ' ' )
                {
                    index++;
                }
            }
            while ( index < text.length() );
            
            offset = offsetTarget + columns[j][i].length();
            
            if ( offsetTarget + lengthMaxes[j] + spacing > w.ws_col - 5 )
            {
                if ( j < columns.size() - 1 )
                {
                    cout << endl;
                }
                
                offset = 0;
            }
            else
            {
                offsetTarget += lengthMaxes[j] + spacing;
            }
        }
        
        cout << endl << endl;
    }
}
