#include "CommandFind.h"
#include "Index.h"

CommandFind::CommandFind()
{
    name = "find";
    description = "Compare query sequences to a reference index";
    argumentString = "index.mash fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
}

int CommandFind::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    Index index;
    
    index.initFromCapnp(arguments[0].c_str());
    
    return 0;
}
