#ifndef INCLUDED_CommandDistance
#define INCLUDED_CommandDistance

#include "Command.h"
#include "Sketch.h"

class CommandDistance : public Command
{
public:
    
    struct CompareInput
    {
        CompareInput(const Sketch & sketchRefNew, Sketch * sketchQueryNew, const std::string fileNew, int kmerSizeNew, float sketchSizeNew, bool concatNew, bool noncanonicalNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            file(fileNew),
            kmerSize(kmerSizeNew),
            sketchSize(sketchSizeNew),
            concat(concatNew),
            noncanonical(noncanonicalNew)
            {}
        
        const Sketch & sketchRef;
        Sketch * sketchQuery;
        std::string nameRef;
        const std::string file;
        int kmerSize;
        int sketchSize;
        bool concat;
        bool noncanonical;
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
float compareSketches(const std::vector<Sketch::hash_t> & hashesSortedRef, const std::vector<Sketch::hash_t> & hashesSortedQuery);

#endif
