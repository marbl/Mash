// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandPaste.h"
#include "Sketch.h"
#include <iostream>
#include "unistd.h"

using std::string;
using std::cerr;
using std::endl;

namespace mash {

CommandPaste::CommandPaste()
: Command()
{
    name = "paste";
    summary = "Create a single sketch file from multiple sketch files.";
    description = "Create a single sketch file from multiple sketch files.";
    argumentString = "<out_prefix> <sketch> [<sketch>] ...";
    
    useOption("help");
    addOption("list", Option(Option::Boolean, "l", "", "Input files are lists of file names.", ""));
}

int CommandPaste::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    bool list = options.at("list").active;
    std::vector<string> files;
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        if ( list )
        {
            splitFile(arguments[i], files);
        }
        else
        {
            files.push_back(arguments[i]);
        }
    }
    
    Sketch sketch;
    std::vector<string> filesGood;
    Sketch::Parameters parameters;
    parameters.parallelism = 1;
    
    for ( int i = 0; i < files.size(); i++ )
    {
        const string & file = files[i];
        
        if ( ! hasSuffix(file, suffixSketch) )
        {
            cerr << "ERROR: The file \"" << file << "\" does not look like a sketch." << endl;
            return 1;
        }
        
        filesGood.push_back(file);
    }
    
	sketch.initFromFiles(filesGood, parameters);
    string out = arguments[0];
    
    if ( ! hasSuffix(out, suffixSketch) )
    {
        out += suffixSketch;
    }

	if( access(out.c_str(), F_OK) != -1 )
	{
		cerr << "ERROR: \"" << out << "\" exists; remove to write." << endl;
		exit(1);
	}
	
    cerr << "Writing " << out << "..." << endl;
    sketch.writeToCapnp(out.c_str());
    
    return 0;
}

} // namespace mash
