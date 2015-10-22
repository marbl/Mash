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
    addOption("header", Option(Option::Boolean, "H", "", "Only show header info. Do not list each sketch.", ""));
}

int CommandInfo::run() const
{
    if ( arguments.size() != 1 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    bool header = options.at("header").active;
    
    const string & file = arguments[0];
    
    if ( ! hasSuffix(file, suffixSketch) )
    {
        cerr << "ERROR: The file \"" << file << "\" does not look like a sketch." << endl;
        return 1;
    }
    
    Sketch sketch;
    
    sketch.initFromCapnp(file.c_str(), header);
    
    cout << "Header:" << endl;
    cout << "  Kmer:                          " << sketch.getKmerSize() << endl;
    cout << "  Target min-hashes per sketch:  " << sketch.getMinHashesPerWindow() << endl;
    cout << "  Canonical kmers:               " << (sketch.getNoncanonical() ? "no" : "yes") << endl;
    
    if ( ! header )
    {
        cout << endl;
        cout << "Sketches (" << sketch.getReferenceCount() << "):" << endl;
        
        vector<vector<string>> columns(4);
        
        columns[0].push_back("Hashes");
        columns[1].push_back("Length");
        columns[2].push_back("ID");
        columns[3].push_back("Comment");
        
        for ( int i = 0; i < sketch.getReferenceCount(); i++ )
        {
            const Sketch::Reference & ref = sketch.getReference(i);
            
            columns[0].push_back(to_string(ref.hashesSorted.size()));
            columns[1].push_back(to_string(ref.length));
            columns[2].push_back(ref.name);
            columns[3].push_back(ref.comment);
        }
        
        printColumns(columns, 2, 2, "-", 0);
    }
    
    return 0;
}
