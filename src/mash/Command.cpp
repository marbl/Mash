#include <iostream>
#include <sys/ioctl.h>

#include "Command.h"

using namespace::std;

Command::Option::Option
(
    Type typeNew,
    std::string identifierNew,
    std::string descriptionNew,
    std::string argumentDefaultNew,
    float argumentMinNew,
    float argumentMaxNew
)
    :
    type(typeNew),
    identifier(identifierNew),
    description(descriptionNew),
    argument(argumentDefaultNew),
    argumentDefault(argumentDefaultNew),
    argumentMin(argumentMinNew),
    argumentMax(argumentMaxNew),
    active(false)
{
}

void Command::Option::setArgument(string argumentNew)
{
    argument = argumentNew;
    
    if ( type == Number || type == Integer )
    {
        bool failed = false;
    
        try
        {
            argumentAsNumber = stof(argument);
        
            if ( argumentMin != argumentMax && (argumentAsNumber < argumentMin || argumentAsNumber > argumentMax) )
            {
                failed = true;
            }
            else if ( type == Integer && int(argumentAsNumber) != argumentAsNumber )
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
            cerr << "ERROR: Argument to -" << identifier << " must be a" << (type == Integer ? "n integer" : " number");
        
            if ( argumentMin != argumentMax )
            {
                cerr << " between " << argumentMin << " and " << argumentMax;
            }
        
            cerr << " (" << argument << " given)" << endl;
            exit(1);
        }
    }
}

void Command::addOption(string name, Option option)
{
    options[name] = option;
    optionNamesByIdentifier[option.identifier] = name;
}

Command::Command()
{
    addAvailableOption("help", Option(Option::Boolean, "h", "Help", ""));
    addAvailableOption("kmer", Option(Option::Integer, "k", "Kmer size. Hashes will be based on strings of this many nucleotides.", "11", 1, 32));
    addAvailableOption("windowed", Option(Option::Boolean, "w", "Windowed", ""));
    addAvailableOption("window", Option(Option::Integer, "l", "Window length. Hashes that are minima in any window of this size will be stored.", "1000"));
    addAvailableOption("mins", Option(Option::Integer, "m", "Min-hashes (not used with -w)", "10000"));
    addAvailableOption("minsWindowed", Option(Option::Integer, "n", "Min-hashes per window (only used with -w).", "10"));
    addAvailableOption("verbose", Option(Option::Boolean, "v", "Verbose", ""));
    addAvailableOption("silent", Option(Option::Boolean, "s", "Silent", ""));
}

void Command::print() const
{
    vector<vector<string>> columns(1);
    
    cout << endl << "Usage:" << endl << endl;
    
    columns[0].push_back("mash " + name + " [options] " + argumentString);
    printColumns(columns);
    
    cout << "Description:" << endl << endl;
    
    columns[0].clear();
    columns[0].push_back(description);
    printColumns(columns);
    
    if ( options.size() == 0 )
    {
        return;
    }
    
    cout << "Options:" << endl << endl;
    
    columns.clear();
    columns.resize(4);
    
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
            case Option::Integer:
                type = "<integer>";
                break;
            case Option::File:
                type = "<path>  ";
                break;
        }
        
        columns[1].push_back(type);
        columns[2].push_back(i->second.argumentDefault);
        columns[3].push_back(i->second.description);
    }
    
    printColumns(columns);
}

int Command::run(int argc, const char ** argv)
{
    for ( int i = 0; i < argc; i++ )
    {
        if ( argv[i][0] == '-' && argv[i][1] != 0 )
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
                
                option.setArgument(argv[i]);
            }
        }
        else
        {
            arguments.push_back(argv[i]);
        }
    }
    
    return run();
}

void Command::useOption(string name)
{
    options[name] = optionsAvailable.at(name);
    optionNamesByIdentifier[options[name].identifier] = name;
}

void Command::addAvailableOption(string name, Option option)
{
    optionsAvailable[name] = option;
}

void printColumns(vector<vector<string>> columns, int indent, int spacing, const char * missing)
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    
    vector<int> lengthMaxes(columns.size(), 0);
    
    for ( int i = 0; i < columns.size(); i++ )
    {
        for ( int j = 0; j < columns[i].size(); j++ )
        {
            if ( columns[i][j].length() == 0 )
            {
                columns[i][j] = missing;
            }
            
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
