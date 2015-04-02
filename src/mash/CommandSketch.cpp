#include "CommandSketch.h"
#include "Index.h"

using namespace::std;

CommandSketch::CommandSketch()
{
    name = "sketch";
    description = "Create a reference index";
    argumentString = "fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
    addOption("kmer", Option(Option::Number, "k", "Kmer size. Hashes will be based on strings of this many nucleotides.", "11"));
    addOption("windowed", Option(Option::Boolean, "w", "Windowed", ""));
    addOption("window", Option(Option::Number, "l", "Window length. Hashes that are minima in any window of this size will be stored.", "1000"));
    addOption("mins", Option(Option::Number, "m", "Min-hashes (not used with -w)", "10000"));
    addOption("minsWindowed", Option(Option::Number, "n", "Min-hashes per window (only used with -w).", "10"));
    addOption("prefix", Option(Option::File, "o", "Output prefix (first input file used if unspecified). The suffix '.mash' will be appended.", ""));
    addOption("verbose", Option(Option::Boolean, "v", "Verbose", ""));
    addOption("silent", Option(Option::Boolean, "s", "Silent", ""));
}

int CommandSketch::run() const
{
    if ( arguments.size() == 0 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int kmer = options.at("kmer").getArgumentAsNumber(1, 32);
    bool windowed = options.at("windowed").active;
    int minHashesPerWindow = windowed ? options.at("minsWindowed").getArgumentAsNumber() : options.at("mins").getArgumentAsNumber();
    int windowSize = options.at("window").getArgumentAsNumber();
    int verbosity = options.at("silent").active ? 0 : options.at("verbose").active ? 2 : 1;
    
    Index index;
    
    index.initFromSequence(arguments, kmer, minHashesPerWindow, windowed, windowSize, false, verbosity);
    
    string prefix = options.at("prefix").argument.length() > 0 ? options.at("prefix").argument : arguments[0];
    
    index.writeToCapnp((prefix + ".mash").c_str());
    
    return 0;
}
