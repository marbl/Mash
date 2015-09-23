#include "CommandPaste.h"
#include "Sketch.h"
#include <iostream>
#include "unistd.h"

using namespace::std;

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
    vector<string> files;
    
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
    
    for ( int i = 0; i < files.size(); i++ )
    {
        const string & file = files[i];
        
        if ( ! hasSuffix(file, suffixSketch) )
        {
            cerr << "ERROR: The file \"" << file << "\" does not look like a sketch." << endl;
            return 1;
        }
        
        if ( i > 0 )
        {
            Sketch sketchTest;
            sketchTest.initFromCapnp(file.c_str(), true);
            
            if ( sketchTest.getKmerSize() != sketch.getKmerSize() )
            {
                cerr << "\nWARNING: The sketch " << file << " has a kmer size (" << sketchTest.getKmerSize() << ") that does not match the first sketch (" << sketch.getKmerSize() << "). This sketch will be skipped.\n\n";
                continue;
            }
            
            if ( sketchTest.getNoncanonical() != sketch.getNoncanonical() )
            {
                cerr << "\nWARNING: The sketch " << file << " is " << (sketchTest.getNoncanonical() ? "noncanonical" : "canonical") << " but the first sketch is not. This sketch will be skipped.\n\n";
                continue;
            }
            
            if ( sketchTest.getMinHashesPerWindow() != sketch.getMinHashesPerWindow() )
            {
                cerr << "\nWARNING: The sketch " << file << " has a target sketch size (" << sketchTest.getMinHashesPerWindow() << ") that does not match the first sketch (" << sketch.getMinHashesPerWindow() << "). This sketch will be skipped.\n\n";
                continue;
            }
        }
        
        sketch.initFromCapnp(file.c_str(), false, i > 0);
    }
    
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
