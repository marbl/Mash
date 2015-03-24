#ifndef INCLUDED_CommandCompare
#define INCLUDED_CommandCompare

#include "Command.h"
#include "Index.h"

class CommandCompare : public Command
{
public:
    
    struct CompareInput
    {
        CompareInput(const Index::Hash_set & minHashesRefNew, const std::string fileNew, int kmerSizeNew, int minsNew)
            :
            minHashesRef(minHashesRefNew),
            file(fileNew),
            kmerSize(kmerSizeNew),
            mins(minsNew)
            {}
        
        const Index::Hash_set & minHashesRef;
        const std::string file;
        int kmerSize;
        int mins;
    };
    
    struct CompareOutput
    {
        float score;
        std::string file;
    };
    
    CommandCompare();
    
    int run() const; // override
    
private:
    
    void writeOutput(CompareOutput * output) const;
};

CommandCompare::CompareOutput * compare(CommandCompare::CompareInput * data);
void getMinHashesForFile(Index::Hash_set & minHashes, const std::string & file, int kmerSize, int mins);

#endif
