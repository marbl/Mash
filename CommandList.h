#ifndef INCLUDED_CommandList
#define INCLUDED_CommandList

#include <map>

#include "Command.h"

class CommandList
{
    std::map<std::string, Command *> commands;
    
public:
    
    ~CommandList();
    
    void addCommand(Command * command);
    void print();
    int run(int argc, const char ** argv);
};

#endif
