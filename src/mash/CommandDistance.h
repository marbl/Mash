#ifndef INCLUDED_CommandDistance
#define INCLUDED_CommandDistance

#include "Command.h"
#include "Index.h"

class CommandDistance : public Command
{
public:
    
    struct CompareInput
    {
        CompareInput(const Index::Hash_set & minHashesRefNew, const std::string nameRefNew, const std::string fileNew, int kmerSizeNew, int minsNew, bool concatNew)
            :
            minHashesRef(minHashesRefNew),
            nameRef(nameRefNew),
            file(fileNew),
            kmerSize(kmerSizeNew),
            mins(minsNew),
            concat(concatNew)
            {}
        
        const Index::Hash_set & minHashesRef;
        std::string nameRef;
        const std::string file;
        int kmerSize;
        int mins;
        bool concat;
    };
    
    struct CompareOutput
    {
        struct PairOutput
        {
            float score;
            std::string file;
        };
        
        std::string nameRef;
        std::vector<PairOutput> pairs;
    };
    
    CommandDistance();
    
    int run() const; // override
    
private:
    
    void writeOutput(CompareOutput * output) const;
};

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * data);
void getMinHashesForFile(Index::Hash_set & minHashes, const std::string & file, int kmerSize, int mins);

#endif
