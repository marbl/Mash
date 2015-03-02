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
		FindData(const Index & indexNew, char * seqIdNew, char * seqNew, uint32_t lengthNew, float thresholdNew)
			:
			index(indexNew),
			length(lengthNew),
			threshold(thresholdNew)
		{
			seq = new char[strlen(seqNew) + 1];
			seqId = new char [strlen(seqIdNew) + 1];
			
			strcpy(seq, seqNew);
			strcpy(seqId, seqIdNew);
		}
		
		~FindData()
		{
			delete [] seq;
			delete [] seqId;
		}
		
		const Index & index;
		char * seqId;
		char * seq;
		uint32_t length;
		float threshold;
	};
	
    CommandFind();
    
    int run() const; // override
};

void * find(void *);

#endif
