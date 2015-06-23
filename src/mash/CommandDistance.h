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
            double score;
            double pValue;
            std::string nameRef;
            std::string nameQuery;
        };
        
        std::vector<PairOutput> pairs;
    };
    
    CommandDistance();
    
    int run() const; // override
    
private:
    
    void writeOutput(CompareOutput * output, bool table, bool log, double pValueMax) const;
};

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * data);
void compareSketches(CommandDistance::CompareOutput::PairOutput & output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, int sketchSize, double kmerSpace);
double pValue(uint32_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint32_t sketchSize);

#endif
