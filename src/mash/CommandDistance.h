#ifndef INCLUDED_CommandDistance
#define INCLUDED_CommandDistance

#include "Command.h"
#include "Sketch.h"

class CommandDistance : public Command
{
public:
    
    struct CompareInput
    {
        CompareInput(const Sketch & sketchRefNew, Sketch * sketchQueryNew, const std::string fileNew, int kmerSizeNew, float errorNew, bool concatNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            file(fileNew),
            kmerSize(kmerSizeNew),
            error(errorNew),
            concat(concatNew)
            {}
        
        const Sketch & sketchRef;
        Sketch * sketchQuery;
        std::string nameRef;
        const std::string file;
        int kmerSize;
        float error;
        bool concat;
    };
    
    struct CompareOutput
    {
        struct PairOutput
        {
            float score;
            std::string nameRef;
            std::string nameQuery;
        };
        
        std::vector<PairOutput> pairs;
    };
    
    CommandDistance();
    
    int run() const; // override
    
private:
    
    void writeOutput(CompareOutput * output) const;
};

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * data);
void getMinHashesForFile(Sketch::Hash_set & minHashes, const std::string & file, int kmerSize, int mins);

#endif
