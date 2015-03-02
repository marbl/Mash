#ifndef INCLUDED_Command
#define INCLUDED_Command

#include <map>
#include <string>
#include <vector>

class Command
{

// A base class for handling command line parsing and usage messages. Specific
// commands can be implemented by subclassing and overriding the run() method.
// The run() method will be called by this base class from run(argc, argv) once
// it has parsed options, which can be accessed from the 'options' map using the
// name strings given when adding them (not 'identifier', which is its parsing
// name to appear after the '-'). Arguments (i.e. operands given after the
// options on the command line) will be available in the 'arguments' vector.

public:
    
    class Option
    {
    public:
    
        enum Type
        {
            Boolean,
            Number,
            File
        }
        type;
        
        std::string identifier;
        std::string description;
        std::string argument;
        std::string argumentDefault;
        
        bool active;
        
        Option() {}
        Option(Type typeNew, std::string identifierNew, std::string descriptionNew, std::string argumentDefaultNew = "");
        
        float getArgumentAsNumber(float min = 0, float max = 0) const;
    };
    
    virtual ~Command() {};
    
    void addOption(std::string name, Option option);
    void print() const;
    virtual int run() const = 0;
    int run(int argc, const char ** argv);
    
    std::string name;
    std::string description;
    std::string argumentString;
    
protected:

    std::map<std::string, Option> options;
    std::vector<std::string> arguments;
    
private:
    
    std::map<std::string, std::string> optionNamesByIdentifier;
};

void printColumns(std::vector<std::vector<std::string>> columns, int indent = 3, int spacing = 3, const char * missing = "-");

#endif
