// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include <iostream>
#include <sys/ioctl.h>
#include <sstream>
#include <fstream>

#include "Command.h"
#include "version.h"

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
    	if ( argument.size() == 0 )
    	{
    		argumentAsNumber = 0;
    		return;
    	}
    	
        bool failed = false;
    	
        try
        {
            argumentAsNumber = stof(argument);
        
            if ( argumentMin != argumentMax && (argumentAsNumber < argumentMin || argumentAsNumber > argumentMax) )
            {
                failed = true;
            }
            else if ( type == Integer && uint64_t(argumentAsNumber) != argumentAsNumber )
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
    else if ( type == Size )
    {
    	if ( argument.size() == 0 )
    	{
    		argumentAsNumber = 0;
    		return;
    	}
    	
    	char suffix = argument[argument.size() - 1];
    	uint64_t factor = 1;
    	
    	if ( suffix < '0' || suffix > '9' )
    	{
    		switch ( suffix )
    		{
    			case 'k':
    			case 'K':
    				factor = 1000;
    				break;
    			case 'm':
    			case 'M':
    				factor = 1000000;
    				break;
    			case 'g':
    			case 'G':
    				factor = 1000000000;
    				break;
    			case 't':
    			case 'T':
    				factor = 1000000000000;
    				break;
    			default:
					cerr << "ERROR: Unrecognized unit (\"" << suffix << "\") in argument to -" << identifier << ". If specified, unit must be one of [kKmMgGtT]." << endl;
					exit(1);
    		}
    		
    		argument.resize(argument.size() - 1);
    	}
    	
    	bool fail = false;
    	
    	try
    	{
    		argumentAsNumber = stof(argument);
    	}
        catch ( const exception & e )
        {
        	fail = true;
        }
        
        if ( argumentAsNumber <= 0 || (uint64_t)argumentAsNumber != argumentAsNumber )
        {
        	fail = true;
        }
        
        if ( fail )
        {
            cerr << "ERROR: Argument to -" << identifier << " must be a whole number, optionally followed by one of [kKmMgGtT]." << endl;
            exit(1);
        }
        
        argumentAsNumber *= factor;
    }
}

void Command::addOption(string name, Option option)
{
    options[name] = option;
    optionNamesByCategory[option.category].push_back(name);
    optionNamesByIdentifier[option.identifier] = name;
}

Command::Command()
{
    addAvailableOption("help", Option(Option::Boolean, "h", "", "Help", ""));
    addAvailableOption("kmer", Option(Option::Integer, "k", "Sketch", "K-mer size. Hashes will be based on strings of this many nucleotides. Canonical nucleotides are used by default (see Alphabet options below).", "21", 1, 32));
    addAvailableOption("windowed", Option(Option::Boolean, "W", "Sketch", "Windowed", ""));
    addAvailableOption("window", Option(Option::Integer, "L", "Window", "Window length. Hashes that are minima in any window of this size will be stored.", "10000"));
    //addAvailableOption("error", Option(Option::Number, "e", "Sketch", "Error bound. The (maximum) number of min-hashes in each sketch will be one divided by this number squared.", "0.05"));
    addAvailableOption("sketchSize", Option(Option::Integer, "s", "Sketch", "Sketch size. Each sketch will have at most this many non-redundant min-hashes.", "1000"));
    addAvailableOption("verbose", Option(Option::Boolean, "v", "Output", "Verbose", ""));
    addAvailableOption("silent", Option(Option::Boolean, "s", "Output", "Silent", ""));
    addAvailableOption("individual", Option(Option::Boolean, "i", "Sketch", "Sketch individual sequences, rather than whole files, e.g. for multi-fastas of single-chromosome genomes or pair-wise gene comparisons.", ""));
    addAvailableOption("warning", Option(Option::Number, "w", "Sketch", "Probability threshold for warning about low k-mer size.", "0.01", 0, 1));
    addAvailableOption("reads", Option(Option::Boolean, "r", "Sketch", "Input is a read set. See Reads options below. Incompatible with -i.", ""));
    addAvailableOption("memory", Option(Option::Size, "b", "Reads", "Use a Bloom filter of this size (raw bytes or with K/M/G/T) to filter out unique k-mers. This is useful if exact filtering with -m uses too much memory. However, some unique k-mers may pass erroneously, and copies cannot be counted beyond 2. Implies -r."));
    addAvailableOption("minCov", Option(Option::Integer, "m", "Reads", "Minimum copies of each k-mer required to pass noise filter for reads. Implies -r.", "1"));
    addAvailableOption("targetCov", Option(Option::Number, "c", "Reads", "Target coverage. Sketching will conclude if this coverage is reached before the end of the input file (estimated by average k-mer multiplicity). Implies -r."));
    addAvailableOption("genome", Option(Option::Size, "g", "Reads", "Genome size. If specified, will be used for p-value calculation instead of an estimated size from k-mer content. Implies -r."));
    addAvailableOption("noncanonical", Option(Option::Boolean, "n", "Alphabet", "Preserve strand (by default, strand is ignored by using canonical DNA k-mers, which are alphabetical minima of forward-reverse pairs). Implied if an alphabet is specified with -a or -z.", ""));
    addAvailableOption("protein", Option(Option::Boolean, "a", "Alphabet", "Use amino acid alphabet (A-Z, except BJOUXZ). Implies -n, -k 9.", ""));
    addAvailableOption("alphabet", Option(Option::String, "z", "Alphabet", "Alphabet to base hashes on (case ignored by default; see -Z). K-mers with other characters will be ignored. Implies -n.", ""));
    addAvailableOption("case", Option(Option::Boolean, "Z", "Alphabet", "Preserve case in k-mers and alphabet (case is ignored by default). Sequence letters whose case is not in the current alphabet will be skipped when sketching.", ""));
    addAvailableOption("threads", Option(Option::Integer, "p", "", "Parallelism. This many threads will be spawned for processing.", "1"));
    addAvailableOption("pacbio", Option(Option::Boolean, "pacbio", "", "Use default settings for PacBio sequences.", ""));
    addAvailableOption("illumina", Option(Option::Boolean, "illumina", "", "Use default settings for Illumina sequences.", ""));
    addAvailableOption("nanopore", Option(Option::Boolean, "nanopore", "", "Use default settings for Oxford Nanopore sequences.", ""));
    addAvailableOption("factor", Option(Option::Number, "f", "Window", "Compression factor", "100"));
    
    addCategory("", "");
    addCategory("Input", "Input");
    addCategory("Output", "Output");
    addCategory("Sketch", "Sketching");
    addCategory("Window", "Sketching (windowed)");
    addCategory("Reads", "Sketching (reads)");
    addCategory("Alphabet", "Sketching (alphabet)");
}

void Command::print() const
{
    vector<vector<string>> columns(1);
    
    cout << endl << "Version: " << version << endl;
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
		if ( *i != "" && optionNamesByCategory.at(*i).size() != 0 )
		{
			dividers.push_back(pair<int, string>(columns[0].size(), "..." + categoryDisplayNames.at(*i) + "..."));
		}
		
		for ( vector<string>::const_iterator j = optionNamesByCategory.at(*i).begin(); j != optionNamesByCategory.at(*i).end(); j++ )
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
					case Option::Size:
						type = "size";
						break;
					case Option::File:
						type = "path";
						break;
					case Option::String:
						type = "text";
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

void Command::useSketchOptions()
{
    useOption("threads");
    useOption("kmer");
    useOption("noncanonical");
    useOption("protein");
    useOption("alphabet");
    useOption("case");
#ifdef COMMAND_FIND
    useOption("windowed");
    useOption("window");
    useOption("factor");
#endif
    useOption("sketchSize");
    useOption("individual");
    useOption("warning");
    useOption("reads");
    useOption("memory");
    useOption("minCov");
    useOption("targetCov");
    useOption("genome");
    //useOption("illumina");
    //useOption("pacbio");
    //useOption("nanopore");
}

void Command::addAvailableOption(string name, Option option)
{
    optionsAvailable[name] = option;
}

void Command::addCategory(string name, string displayName)
{
	// keep track of category order
	//
	if ( categoryDisplayNames.count(name) == 0 )
	{
		categories.push_back(name);
	    categoryDisplayNames[name] = displayName;
	    optionNamesByCategory[name] = vector<string>();
	}
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
