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
            Integer,
            File
        }
        type;
        
        std::string identifier;
        std::string description;
        std::string argument;
        std::string argumentDefault;
        
        float argumentAsNumber;
        float argumentMin;
        float argumentMax;
        
        bool active;
        bool changed;
        
        Option() {}
        Option(Type typeNew, std::string identifierNew, std::string descriptionNew, std::string argumentDefaultNew = "", float argumentMinNew = 0, float argumentMaxNew = 0);
        
        float getArgumentAsNumber() const {return argumentAsNumber;}
        void setArgument(std::string argumentNew);
    };
    
    Command();
    virtual ~Command() {};
    
    void addOption(std::string name, Option option);
    void print() const;
    virtual int run() const = 0;
    int run(int argc, const char ** argv);
    void useOption(std::string name);
    
    std::string name;
    std::string description;
    std::string argumentString;
    
protected:

    std::map<std::string, Option> options;
    std::map<std::string, Option> optionsAvailable;
    std::vector<std::string> arguments;
    
private:
    
    void addAvailableOption(std::string name, Option option);
    
    std::map<std::string, std::string> optionNamesByIdentifier;
};

void splitFile(const std::string & file, std::vector<std::string> & lines);
void printColumns(std::vector<std::vector<std::string>> columns, int indent = 3, int spacing = 3, const char * missing = "-");

#endif
