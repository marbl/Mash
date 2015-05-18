#include "CommandSketch.h"
#include "Sketch.h"
#include <iostream>

using namespace::std;

CommandSketch::CommandSketch()
: Command()
{
    name = "sketch";
    description = "Create a reference sketch";
    argumentString = "fast(a|q)[.gz] ...";
    
    useOption("help");
    useOption("kmer");
    useOption("windowed");
    useOption("window");
    useOption("factor");
    useOption("verbose");
    useOption("silent");
    useOption("concat");
    useOption("illumina");
    useOption("pacbio");
    useOption("nanopore");
    addOption("prefix", Option(Option::File, "o", "Output prefix (first input file used if unspecified). The suffix '.msh' will be appended.", ""));
}

int CommandSketch::run() const
{
    if ( arguments.size() == 0 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int kmer = options.at("kmer").getArgumentAsNumber();
    bool windowed = options.at("windowed").active;
    float factor = options.at("factor").getArgumentAsNumber();
    int windowSize = options.at("window").getArgumentAsNumber();
    int verbosity = options.at("silent").active ? 0 : options.at("verbose").active ? 2 : 1;
    bool concat = options.at("concat").active;
    
    if ( concat && windowed )
    {
        cerr << "ERROR: " << options.at("concat").identifier << " and " << options.at("windowed").identifier << " are incompatible." << endl;
        return 1;
    }
    
    Sketch sketch;
    
    sketch.initFromSequence(arguments, kmer, factor, windowed, windowSize, concat, verbosity);
    
    string prefix = options.at("prefix").argument.length() > 0 ? options.at("prefix").argument : arguments[0];
    
    //sketch.writeToCapnp((prefix + (windowed ? suffixSketchWindowed : suffixSketch)).c_str());
    
    return 0;
}
