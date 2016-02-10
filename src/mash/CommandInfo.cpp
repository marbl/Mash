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
    addOption("header", Option(Option::Boolean, "H", "", "Only show header info. Do not list each sketch. Incompatible with -t", ""));
    addOption("tabular", Option(Option::Boolean, "t", "", "Tabular output (rather than padded), with no header. Incompatible with -H.", ""));
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
    
    if ( header && tabular )
    {
    	cerr << "ERROR: The options -H and -t are incompatible." << endl;
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
    
    if ( tabular )
    {
    	cout << "#Hashes\tLength\tID\tComment" << endl;
    }
    else
    {
		cout << "Header:" << endl;
		cout << "  Kmer:                          " << sketch.getKmerSize() << endl;
		cout << "  Target min-hashes per sketch:  " << sketch.getMinHashesPerWindow() << endl;
		cout << "  Canonical kmers:               " << (sketch.getNoncanonical() ? "no" : "yes") << endl;
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
