// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandList
#define INCLUDED_CommandList

#include <map>

#include "Command.h"

namespace mash {

class CommandList
{
    std::map<std::string, Command *> commands;
    
public:
    
    CommandList(std::string nameNew);
    ~CommandList();
    
    void addCommand(Command * command);
    void print();
    int run(int argc, const char ** argv);
    
private:
    
    void showLicense();
    
    std::string name;
};

} // namespace mash

#endif
