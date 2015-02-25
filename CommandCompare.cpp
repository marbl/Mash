#include "CommandCompare.h"
#include "Index.h"
#include <iostream>

using namespace::std;

CommandCompare::CommandCompare()
{
    name = "compare";
    description = "Compute the global distance of each input sequence to the reference index. The score is computed as 2*c/(m+n), where m is the min-hash count for the reference, n is the min-hash count for the query, and c is the number of min-hashes in common.";
    argumentString = "reference.mash fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
}

int CommandCompare::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    Index indexRef;
    indexRef.initFromCapnp(arguments[0].c_str());
    
    int l;
    int count = 0;
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        Index indexQuery;
        indexQuery.initFromSequence(vector<string>(1, arguments[i]), indexRef.getKmerSize(), indexRef.getCompressionFactor());
        
        cout << compare(indexRef, indexQuery) << '\t' << arguments[i] << endl;
    }
    
    return 0;
}

float compare(const Index & indexRef, const Index & indexQuery)
{
    int common = 0;
    
    for ( Index::LociByHash_umap::const_iterator i = indexQuery.getLociByHash().begin(); i != indexQuery.getLociByHash().end(); i++ )
    {
        if ( indexRef.getLociByHash().count(i->first) != 0 )
        {
            common++;
        }
    }
    
    return float(common) * 2 / (indexRef.getLociByHash().size() + indexQuery.getLociByHash().size());
}
