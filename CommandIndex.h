#ifndef INCLUDED_CommandIndex
#define INCLUDED_CommandIndex

#include "Command.h"

class CommandIndex : public Command
{
public:

	CommandIndex(std::string nameNew, std::string descriptionNew, std::string argumentStringNew)
		:
		Command(nameNew, descriptionNew, argumentStringNew)
		{}
	
	int run() const;
};

#endif
