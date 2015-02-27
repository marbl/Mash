#ifndef INCLUDED_CommandFind
#define INCLUDED_CommandFind

#include "Command.h"
#include "Index.h"
#include <string.h>

class CommandFind : public Command
{
public:
	
	struct FindData
	{
		FindData(const Index & indexNew, char * seqNew, uint32_t lengthNew, float thresholdNew)
			:
			index(indexNew),
			length(lengthNew),
			threshold(thresholdNew)
		{
			seq = new char[strlen(seqNew) + 1];
			
			strcpy(seq, seqNew);
		}
		
		~FindData()
		{
			delete [] seq;
		}
		
		const Index & index;
		char * seq;
		uint32_t length;
		float threshold;
	};
	
    CommandFind();
    
    int run() const; // override
};

void * find(void *);

#endif
