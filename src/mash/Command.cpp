#include <iostream>
#include <sys/ioctl.h>
#include <sstream>
#include <fstream>

#include "Command.h"

using namespace::std;

Command::Option::Option
(
    Type typeNew,
    string identifierNew,
    string categoryNew,
    string descriptionNew,
    string argumentDefaultNew,
    float argumentMinNew,
    float argumentMaxNew
)
    :
    type(typeNew),
    identifier(identifierNew),
    category(categoryNew),
    description(descriptionNew),
    argumentDefault(argumentDefaultNew),
    argumentMin(argumentMinNew),
    argumentMax(argumentMaxNew),
    active(false)
{
    setArgument(argumentDefault);
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
	// keep track of category order
	//
	if ( optionNamesByCategory.count(option.category) == 0 )
	{
		categories.push_back(option.category);
	}
	
    options[name] = option;
    optionNamesByCategory[option.category].insert(name);
    optionNamesByIdentifier[option.identifier] = name;
}

Command::Command()
{
    addAvailableOption("help", Option(Option::Boolean, "h", "", "Help", ""));
    addAvailableOption("kmer", Option(Option::Integer, "k", "Sketch", "Kmer size. Hashes will be based on strings of this many nucleotides.", "15", 1, 32));
    addAvailableOption("windowed", Option(Option::Boolean, "w", "Sketch", "Windowed", ""));
    addAvailableOption("window", Option(Option::Integer, "l", "Sketch", "Window length. Hashes that are minima in any window of this size will be stored.", "1000"));
    addAvailableOption("error", Option(Option::Number, "e", "Sketch", "Error bound. The (maximum) number of min-hashes in each sketch will be one divided by this number squared.", "0.05"));
    addAvailableOption("sketchSize", Option(Option::Integer, "s", "Sketch", "Sketch size. Each sketch will have at most this many non-redundant min-hashes.", "10000"));
    addAvailableOption("verbose", Option(Option::Boolean, "v", "Output", "Verbose", ""));
    addAvailableOption("silent", Option(Option::Boolean, "s", "Output", "Silent", ""));
    addAvailableOption("individual", Option(Option::Boolean, "i", "Sketch", "Sketch individual sequences, rather than whole files.", ""));
    addAvailableOption("unique", Option(Option::Boolean, "u", "Sketch", "Remove (most) unique kmers using a Bloom Filter. This is useful for reducing noise from sequencing errors in read sets. See Bloom filter options below. Implies -f.", ""));
    addAvailableOption("genome", Option(Option::Integer, "g", "Bloom", "Expected genome size (Mb). Helps pick the Bloom Filter size. Should be within an order of magnitude of the true size. Implies -u.", "5"));
    addAvailableOption("memory", Option(Option::Integer, "m", "Bloom", "Maximum Bloom Filter memory usage (GB). More memory will allow more thorough detection of unique kmers, so this should be as high as is practical for the computing environment (though it may not actually be used). Implies -u.", "1", 1, 1024));
    addAvailableOption("bloomError", Option(Option::Number, "e", "Bloom", "Target false-negative rate for filtering unique kmers with -u.", "0.1", 0, 1));
    addAvailableOption("noncanonical", Option(Option::Boolean, "n", "Sketch", "Non-canonical. By default, canonical DNA kmers (alphabetical minima of forward-reverse pairs) are used, and kmers with non-acgtACGT characters are ignored. This option uses kmers as they appear and allows all characters.", ""));
    addAvailableOption("threads", Option(Option::Integer, "p", "", "Parallelism. This many threads will be spawned, each one handling one query sequence at a time.", "1"));
    addAvailableOption("pacbio", Option(Option::Boolean, "pacbio", "", "Use default settings for PacBio sequences.", ""));
    addAvailableOption("illumina", Option(Option::Boolean, "illumina", "", "Use default settings for Illumina sequences.", ""));
    addAvailableOption("nanopore", Option(Option::Boolean, "nanopore", "", "Use default settings for Oxford Nanopore sequences.", ""));
    
    categoryDisplayNames["Input"] = "Input";
    categoryDisplayNames["Output"] = "Output";
    categoryDisplayNames["Sketch"] = "Sketching";
    categoryDisplayNames["Bloom"] = "Sketching (Bloom filter)";
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
	columns.resize(2);

	columns[0].push_back("Option");
	columns[1].push_back("Description (range) [default]");
    
    vector<pair<int, string> > dividers;
    
    for ( vector<string>::const_iterator i = categories.begin(); i != categories.end(); i++ )
    {
		if ( *i != "" )
		{
			dividers.push_back(pair<int, string>(columns[0].size(), "..." + categoryDisplayNames.at(*i) + "..."));
		}
		
		for ( set<string>::const_iterator j = optionNamesByCategory.at(*i).begin(); j != optionNamesByCategory.at(*i).end(); j++ )
		{
			const Option & option = options.at(*j);
		
			string optionString = "-" + option.identifier;
		
			string range;
		
			if ( option.type != Option::Boolean )
			{
				string type;
			
				switch ( option.type )
				{
					case Option::Boolean:
						break;
					case Option::Number:
						type = "num";
						break;
					case Option::Integer:
						type = "int";
						break;
					case Option::File:
						type = "path";
						break;
				}
			
				optionString += " <" + type + ">";
			}
		
			columns[0].push_back(optionString);
		
			string descString = option.description;
		
			if ( option.argumentMin != option.argumentMax )
			{
				stringstream stringMin;
				stringstream stringMax;
			
				if ( option.type == Option::Integer )
				{
					stringMin << int(option.argumentMin);
					stringMax << int(option.argumentMax);
				}
				else
				{
					stringMin << option.argumentMin;
					stringMax << option.argumentMax;
				}
			
				descString += " (" + stringMin.str() + "-" + stringMax.str() + ")";
			}
		
			if ( option.argumentDefault != "" )
			{
				descString += " [" + option.argumentDefault + "]";
			}
		
			columns[1].push_back(descString);
		}
	}
	
	printColumns(columns, dividers);
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
            
            option.active = true;
            
            if ( option.type != Option::Boolean )
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
    addOption(name, optionsAvailable.at(name));
}

void Command::addAvailableOption(string name, Option option)
{
    optionsAvailable[name] = option;
}

void splitFile(const string & file, vector<string> & lines)
{
    string line;
    
    ifstream in(file);
    
    if ( in.fail() )
    {
    	cerr << "ERROR: Could not open " << file << ".\n";
    	exit(1);
    }
    
    while ( getline(in, line) )
    {
        lines.push_back(line);
    }
}

void printColumns(const vector<vector<string>> & columns, int indent, int spacing, const char * missing, int max)
{
	printColumns(columns, vector<pair<int, string>>(), indent, spacing, missing, max);
}

void printColumns(const vector<vector<string>> & columns, const vector<pair<int, string>> & dividers, int indent, int spacing, const char * missing, int max)
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    
    int cols = w.ws_col;
    
    if ( max != 0 && max < cols )
    {
    	cols = max;
    }
    
    int dividerCurrent = 0;
    
    vector<int> lengthMaxes(columns.size(), 0);
    
    for ( int i = 0; i < columns.size(); i++ )
    {
        for ( int j = 0; j < columns[i].size(); j++ )
        {
            int length = columns[i][j].length();
            
            if ( length == 0 )
            {
                length = 1;
            }
            
            if ( length > lengthMaxes[i] )
            {
                lengthMaxes[i] = length;
            }
        }
    }
    
    for ( int i = 0; i < columns[0].size(); i++ )
    {
        int offset = 0;
        int offsetTarget = indent;
        
        if ( dividerCurrent < dividers.size() && i == dividers.at(dividerCurrent).first )
        {
        	cout << dividers.at(dividerCurrent).second << endl << endl;
        	dividerCurrent++;
        }
        
        for ( int j = 0; j < columns.size(); j++ )
        {
            for ( int k = offset; k < offsetTarget; k++ )
            {
                cout << ' ';
            }
            
            int index = 0;
            string text = columns[j][i];
            
            if ( text == "" )
            {
            	text = missing;
            }
            
            do
            {
                int length = text.length() - index;
                
                if ( length + offsetTarget > cols )
                {
                    length = cols - offsetTarget;
                
                    while ( text[index + length] != ' ' && length > 0 )
                    {
                        length--;
                    }
                }
                
                if ( length == 0 )
                {
                    length = cols - offsetTarget;
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
            
            if ( offsetTarget + lengthMaxes[j] + spacing > cols - 5 )
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
