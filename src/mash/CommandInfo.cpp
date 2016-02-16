// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandInfo.h"
#include "Sketch.h"
#include <iostream>

using namespace::std;

CommandInfo::CommandInfo()
: Command()
{
    name = "info";
    summary = "Display information about sketch files.";
    description = "Display information about sketch files.";
    argumentString = "<sketch>";
    
    useOption("help");
    addOption("header", Option(Option::Boolean, "H", "", "Only show header info. Do not list each sketch. Incompatible with -t and -c.", ""));
    addOption("tabular", Option(Option::Boolean, "t", "", "Tabular output (rather than padded), with no header. Incompatible with -H and -c.", ""));
    addOption("counts", Option(Option::Boolean, "c", "", "Show hash count histograms for each sketch. Incompatible with -H and -t.", ""));
}

int CommandInfo::run() const
{
    if ( arguments.size() != 1 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    bool header = options.at("header").active;
    bool tabular = options.at("tabular").active;
    bool counts = options.at("counts").active;
    
    if ( header && tabular )
    {
    	cerr << "ERROR: The options -H and -t are incompatible." << endl;
    	return 1;
    }
    
    if ( header && counts )
    {
    	cerr << "ERROR: The options -H and -c are incompatible." << endl;
    	return 1;
    }
    
    if ( tabular && counts )
    {
    	cerr << "ERROR: The options -t and -c are incompatible." << endl;
    	return 1;
    }
    
    const string & file = arguments[0];
    
    if ( ! hasSuffix(file, suffixSketch) )
    {
        cerr << "ERROR: The file \"" << file << "\" does not look like a sketch." << endl;
        return 1;
    }
    
    Sketch sketch;
    Sketch::Parameters params;
    params.parallelism = 1;
    
    sketch.initFromFiles(arguments, params);
    
    if ( counts )
    {
    	return printCounts(sketch);
    }
    
    if ( tabular )
    {
    	cout << "#Hashes\tLength\tID\tComment" << endl;
    }
    else
    {
    	string alphabet;
    	sketch.getAlphabetAsString(alphabet);
    	
		cout << "Header:" << endl;
		cout << "  K-mer size:                    " << sketch.getKmerSize() << " (" << (sketch.getUse64() ? "64" : "32") << "-bit hashes)" << endl;
		cout << "  Alphabet:                      " << alphabet << (sketch.getNoncanonical() ? "" : " (canonical)") << (sketch.getPreserveCase() ? " (case-sensitive)" : "") << endl;
		cout << "  Target min-hashes per sketch:  " << sketch.getMinHashesPerWindow() << endl;
		cout << "  Sketches:                      " << sketch.getReferenceCount() << endl;
	}
	
    if ( ! header )
    {
        vector<vector<string>> columns(4);
        
        if ( ! tabular )
        {
			cout << endl;
			cout << "Sketches:" << endl;
		
			columns[0].push_back("Hashes");
			columns[1].push_back("Length");
			columns[2].push_back("ID");
			columns[3].push_back("Comment");
		}
        
        for ( uint64_t i = 0; i < sketch.getReferenceCount(); i++ )
        {
            const Sketch::Reference & ref = sketch.getReference(i);
            
            if ( tabular )
            {
            	cout
            		<< ref.hashesSorted.size() << '\t'
            		<< ref.length << '\t'
            		<< ref.name << '\t'
            		<< ref.comment << endl;
            }
            else
            {
				columns[0].push_back(to_string(ref.hashesSorted.size()));
				columns[1].push_back(to_string(ref.length));
				columns[2].push_back(ref.name);
				columns[3].push_back(ref.comment);
			}
        }
        
        if ( ! tabular )
        {
	        printColumns(columns, 2, 2, "-", 0);
	    }
    }
    
    return 0;
}

int CommandInfo::printCounts(const Sketch & sketch) const
{
	if ( sketch.getReferenceCount() == 0 )
	{
		cerr << "ERROR: Sketch file contains no sketches" << endl;
		return 1;
	}
	
	if ( ! sketch.hasHashCounts() )
	{
		cerr << "ERROR: Sketch file does not have hash counts. Re-sketch with counting enabled to use this feature." << endl;
		return 1;
	}
	
	cout << "#Sketch\tBin\tFrequency" << endl;
	
	map<uint32_t, uint64_t> histogram;
	
	for ( uint64_t i = 0; i < sketch.getReferenceCount(); i++ )
	{
		const string & name = sketch.getReference(i).name;
		
		sketch.getReferenceHistogram(i, histogram);
		
		for ( map<uint32_t, uint64_t>::const_iterator j = histogram.begin(); j != histogram.end(); j++ )
		{
			cout << name << '\t' << j->first << '\t' << j->second << endl;
		}
	}
	
	return 0;
}
