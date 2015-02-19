#include "CommandIndex.h"
#include "Index.h"

using namespace::std;

CommandIndex::CommandIndex()
{
	name = "index";
	description = "Create a reference index";
	argumentString = "fast(a|q)[.gz] ...";
	
	addOption("kmer", Option(Option::Number, "k", "Kmer size", "11"));
	addOption("factor", Option(Option::Number, "c", "Compression factor", "100"));
	addOption("hashes", Option(Option::Number, "m", "Number of hashes", "10"));
	addOption("prefix", Option(Option::File, "p", "Output prefix (first input file used if unspecified). The suffix'.mash' will be appended.", ""));
}

int CommandIndex::run() const
{
	if ( arguments.size() == 0 )
	{
		print();
		return 0;
	}
	
	int kmer = options.at("kmer").getArgumentAsNumber();
	int factor = options.at("factor").getArgumentAsNumber();
	int mins = options.at("hashes").getArgumentAsNumber();
	
	Index index;
	
	index.initFromSequence(arguments[0].c_str(), kmer, mins, factor);
	index.writeToCapnp((arguments[0] + ".mash").c_str());
	
	return 0;
}
