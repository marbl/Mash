#ifndef INCLUDED_CommandIndex
#define INCLUDED_CommandIndex

#include "Command.h"

class CommandIndex : public Command
{
public:

	CommandIndex(std::string nameNew, std::string descriptionNew)
		:
		Command(nameNew, descriptionNew)
		{}
	
	int run() const;
};

#endif
