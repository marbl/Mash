#include "CommandPaste.h"
#include "Sketch.h"
#include <iostream>

using namespace::std;

CommandPaste::CommandPaste()
: Command()
{
    name = "paste";
    description = "Create a single sketch file from multiple sketch files.";
    argumentString = "<out_prefix> <sketch> [<sketch>] ...";
    
    useOption("help");
}

int CommandPaste::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    Sketch sketch;
    
    sketch.initFromCapnp(arguments[1].c_str());
    
    for ( int i = 2; i < arguments.size(); i++ )
    {
        const string & file = arguments[i];
        
        if ( ! hasSuffix(file, suffixSketch) )
        {
            cerr << "ERROR: The file \"" << file << "\" does not look like a sketch." << endl;
            return 1;
        }
        
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
        
        sketch.initFromCapnp(file.c_str(), false, true);
    }
    
    string out = arguments[0] + suffixSketch;
    
    cerr << "Writing " << out << "..." << endl;
    sketch.writeToCapnp(out.c_str());
    
    return 0;
}
