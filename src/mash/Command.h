// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_Command
#define INCLUDED_Command

#include <map>
#include <string>
#include <vector>
#include <set>

namespace mash {

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
            Size,
            File,
            String
        }
        type;
        
        std::string category;
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
        Option(Type typeNew, std::string identifierNew, std::string categoryNew, std::string descriptionNew, std::string argumentDefaultNew = "", float argumentMinNew = 0, float argumentMaxNew = 0);
        
        float getArgumentAsNumber() const {return argumentAsNumber;}
        void setArgument(std::string argumentNew);
    };
    
    Command();
    virtual ~Command() {};
    
    void addOption(std::string name, Option option);
    const Option & getOption(std::string name) const;
    inline bool hasOption(std::string name) const {return options.count(name);}
    void print() const;
    virtual int run() const = 0;
    int run(int argc, const char ** argv);
    
    std::string name;
    std::string summary;
    std::string description;
    std::string argumentString;
    
protected:
	
    void useOption(std::string name);
	void useSketchOptions();
	
    std::map<std::string, Option> options;
    std::map<std::string, Option> optionsAvailable;
    std::vector<std::string> arguments;
    
private:
    
    void addAvailableOption(std::string name, Option option);
	void addCategory(std::string name, std::string displayName);
    
    std::map<std::string, std::string> optionNamesByIdentifier;
    std::map<std::string, std::vector<std::string> > optionNamesByCategory;
    std::vector<std::string> categories;
    std::map<std::string, std::string> categoryDisplayNames;
};

inline const Command::Option & Command::getOption(std::string name) const {return options.at(name);}
void splitFile(const std::string & file, std::vector<std::string> & lines);
void printColumns(const std::vector<std::vector<std::string>> & columns, int indent = 2, int spacing = 2, const char * missing = "-", int max = 80);
void printColumns(const std::vector<std::vector<std::string>> & columns, const std::vector<std::pair<int, std::string>> & dividers, int indent = 2, int spacing = 2, const char * missing = "-", int max = 80);

} // namespace mash

#endif
