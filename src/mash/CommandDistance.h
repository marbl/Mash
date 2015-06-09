#ifndef INCLUDED_CommandDistance
#define INCLUDED_CommandDistance

#include "Command.h"
#include "Sketch.h"

class CommandDistance : public Command
{
public:
    
    struct CompareInput
    {
        CompareInput(const Sketch & sketchRefNew, Sketch * sketchQueryNew, const std::string fileNew, const Sketch::Parameters & parametersNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            file(fileNew),
            parameters(parametersNew)
            {}
        
        const Sketch & sketchRef;
        Sketch * sketchQuery;
        std::string nameRef;
        const std::string file;
        const Sketch::Parameters & parameters;
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
float compareSketches(const HashList & hashesSortedRef, const HashList & hashesSortedQuery, int targetSizeRef, int targetSizeQuery);

#endif
