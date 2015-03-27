#include "CommandIndex.h"
#include "Index.h"

using namespace::std;

CommandIndex::CommandIndex()
{
    name = "index";
    description = "Create a reference index";
    argumentString = "fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
    addOption("kmer", Option(Option::Number, "k", "Kmer size. Hashes will be based on strings of this many nucleotides.", "11"));
    addOption("window", Option(Option::Number, "w", "Window size. Hashes that are minima in any window of this size will be stored.", "1000"));
    addOption("mins", Option(Option::Number, "m", "Min-hashes per window.", "10"));
    addOption("prefix", Option(Option::File, "o", "Output prefix (first input file used if unspecified). The suffix '.mash' will be appended.", ""));
    addOption("verbose", Option(Option::Boolean, "v", "Verbose", ""));
}

int CommandIndex::run() const
{
    if ( arguments.size() == 0 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int kmer = options.at("kmer").getArgumentAsNumber(1, 32);
    int minHashesPerWindow = options.at("mins").getArgumentAsNumber();
    int windowSize = options.at("window").getArgumentAsNumber();
    int verbosity = options.at("verbose").active ? 2 : 1;
    
    Index index;
    
    index.initFromSequence(arguments, kmer, minHashesPerWindow, windowSize, verbosity);
    
    string prefix = options.at("prefix").argument.length() > 0 ? options.at("prefix").argument : arguments[0];
    
    index.writeToCapnp((prefix + ".mash").c_str());
    
    return 0;
}
