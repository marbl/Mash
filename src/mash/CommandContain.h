#ifndef INCLUDED_CommandContain
#define INCLUDED_CommandContain

#include "Command.h"
#include "Sketch.h"

class CommandContain : public Command
{
public:
    
    struct ContainInput
    {
        ContainInput(const Sketch & sketchRefNew, Sketch * sketchQueryNew, const std::string fileNew, int kmerSizeNew, float sketchSizeNew, bool concatNew, bool noncanonicalNew)
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
    
    struct ContainOutput
    {
        struct PairOutput
        {
            float score;
            std::string nameRef;
            std::string nameQuery;
        };
        
        std::vector<PairOutput> pairs;
    };
    
    CommandContain();
    
    int run() const; // override
    
private:
    
    void writeOutput(ContainOutput * output) const;
};

CommandContain::ContainOutput * contain(CommandContain::ContainInput * data);
float containSketches(const std::vector<Sketch::hash_t> & hashesSortedRef, const std::vector<Sketch::hash_t> & hashesSortedQuery);

#endif
