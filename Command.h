#ifndef INCLUDED_Command
#define INCLUDED_Command

#include <map>
#include <string>
#include <vector>

class Command
{
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
		
		bool active;
		
		Option() {}
		Option(Type typeNew, std::string identifierNew, std::string descriptionNew, std::string argumentDefault = "")
			:
			type(typeNew),
			identifier(identifierNew),
			description(descriptionNew),
			argument(argumentDefault),
			active(false)
			{}
		
		float getArgumentAsNumber(float min = 0, float max = 0) const;
	};
	
	Command(std::string nameNew, std::string descriptionNew)
		:
		name(nameNew),
		description(descriptionNew)
		{}
	virtual ~Command() {};
	
	void addOption(std::string name, Option option);
	void print() const;
	virtual int run() const = 0;
	int run(int argc, const char ** argv);
	
	std::string name;
	std::string description;
	
protected:

	std::map<std::string, Option> options;
	std::vector<std::string> arguments;
	
private:
	
	std::map<std::string, std::string> optionNamesByIdentifier;
};

#endif
