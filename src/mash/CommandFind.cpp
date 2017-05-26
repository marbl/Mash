// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandFind.h"
#include "Sketch.h"
#include <zlib.h>
#include "kseq.h"
#include <iostream>
#include <set>
#include <unordered_set>
#include "ThreadPool.h"
#include "sketchParameterSetup.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace mash {

KSEQ_INIT(gzFile, gzread)

CommandFind::CommandFind()
: Command()
{
    name = "find";
    summary = "Find regions of references that have similarity to query sequences.";
    description = "Compare query sequences to a reference. <reference> can be a fasta file (gzipped or not) or a mash windowed sketch file (.msw). If it is fasta and a sketch file does not exist (or is out of date), one will be written with the current options and used for future runs if possible. <query> can be fasta or fastq, gzipped or not. Multiple query files can be provided, or \"-\" can be given to read from standard input.";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    addOption("threshold", Option(Option::Number, "t", "Output", "Threshold. This fraction of the query sequence's min-hashes must appear in a query-sized window of a reference sequence for the match to be reported.", "0.2", 0.0, 1.0));
    addOption("best", Option(Option::Integer, "b", "Output", "Best hit count. This many of the best hits will be reported (0 to report all hits). Score ties are broken by keeping the hit to the earlier reference or to the left-most position.", "0"));
    addOption("self", Option(Option::Boolean, "self", "Output", "Ignore self matches if query ID appears in reference.", ""));
    useSketchOptions();
}

int CommandFind::run() const
{
    if ( arguments.size() < 2 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    float threshold = options.at("threshold").getArgumentAsNumber();
    int threads = options.at("threads").getArgumentAsNumber();
    int best = options.at("best").getArgumentAsNumber();
    bool selfMatches = ! options.at("self").active;
    
	Sketch::Parameters params;
	
    if ( sketchParameterSetup(params, *(Command *)this) )
    {
    	return 1;
    }
    
    params.windowed = true;
    Sketch sketch;
    const string & fileReference = arguments[0];
    
    if ( hasSuffix(fileReference, suffixSketch ) )
    {
        cerr << "ERROR: Reference (" << fileReference << ") looks like a sketch but is not windowed.\n";
        return 1;
    }
    
    if ( hasSuffix(fileReference, suffixSketchWindowed) )
    {
        if ( options.at("kmer").active || options.at("sketchSize").active || options.at("window").active )
        {
            cerr << "ERROR: The options " << options.at("kmer").identifier << ", " << options.at("sketchSize").identifier << " and " << options.at("window").identifier << " cannot be used when a sketch is provided; these are inherited from the sketch.\n";
            return 1;
        }
    }
    else
    {
        int kmerSize = options.at("kmer").getArgumentAsNumber();
        float factor = options.at("factor").getArgumentAsNumber();
        int windowSize = options.at("window").getArgumentAsNumber();
        int mins = windowSize / factor;
        
		//cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;
		cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...\n";
		
		params.kmerSize = kmerSize;
		params.minHashesPerWindow = mins;
		params.windowed = true;
		params.windowSize = windowSize;
		params.parallelism = threads;
		
		/*
		if ( sketch.writeToFile() )
		{
			cerr << "Sketch saved for subsequent runs." << endl;
		}
		else
		{
			cerr << "The sketch for " << fileReference << " could not be saved; it will be sketched again next time." << endl;
		}*/
    }
    
	vector<string> refArgVector;
	refArgVector.push_back(fileReference);
	
	sketch.initFromFiles(refArgVector, params);
	
    int l;
    int count = 0;
    
    ThreadPool<FindInput, FindOutput> threadPool(find, threads);
    
    for ( int i = 1; i < arguments.size(); i++ )
    {
        FILE * inStream = 0;
        
        if ( arguments[i] == "-" )
        {
            inStream = stdin;
        }
        else
        {
            inStream = fopen(arguments[i].c_str(), "r");
            
            if ( inStream == NULL )
            {
                cerr << "ERROR: could not open " << arguments[i] << " for reading." << endl;
                exit(1);
            }
        }
        
        gzFile fp = gzdopen(fileno(inStream), "r");
        kseq_t *seq = kseq_init(fp);
        
        while ((l = kseq_read(seq)) >= 0)
        {
            if ( l < sketch.getKmerSize() )
            {
                continue;
            }
            
            //printf("Query name: %s\tlength: %d\n\n", seq->name.s, l);
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            
            threadPool.runWhenThreadAvailable(new FindInput(sketch, seq->name.s, seq->seq.s, l, threshold, best, selfMatches));
            
            while ( threadPool.outputAvailable() )
            {
                //cout << "popping\n";
                writeOutput(sketch, threadPool.popOutputWhenAvailable());
            }
        }
        
        if ( l != -1 )
        {
            printf("ERROR: return value: %d\n", l);
            return 1;
        }
        
        kseq_destroy(seq);
        gzclose(fp);
        fclose(inStream);
    }
    
    while ( threadPool.running() )
    {
        writeOutput(sketch, threadPool.popOutputWhenAvailable());
    }
    
    return 0;
}

void CommandFind::writeOutput(const Sketch & sketch, FindOutput * output) const
{
    //cout << output->seqId << endl;
    
    // reverse the order befor printing
    //
    vector<FindOutput::Hit> hits;
    //
    while ( output->hits.size() > 0 )
    {
        hits.push_back(output->hits.top());
        output->hits.pop();
    }
    
    for ( int i = hits.size() - 1; i >= 0; i-- )
    {
        const FindOutput::Hit & hit = hits.at(i);
        cout <<
            output->seqId << '\t' <<
            sketch.getReference(hit.ref).name << '\t' <<
            hit.start << '\t' <<
            hit.end << '\t' <<
            (hit.minusStrand ? '-' : '+') << '\t' <<
            hit.score << endl;
    }
    
    delete output;
}

CommandFind::FindOutput * find(CommandFind::FindInput * data)
{
    CommandFind::FindOutput * output = new CommandFind::FindOutput();
    
    // uppercase the entire sequence in place
    //
    char * seq = data->seq;
    //
    for ( int i = 0; i < data->length; i++ )
    {
        if ( seq[i] > 90 )
        {
            seq[i] -= 32;
        }
    }
    
    findPerStrand(data, output, false);
    findPerStrand(data, output, true);
    
    return output;
}

void findPerStrand(const CommandFind::FindInput * input, CommandFind::FindOutput * output, bool minusStrand)
{
    typedef std::unordered_map < uint32_t, std::set<uint32_t> > PositionsBySequence_umap;
    
    bool verbose = false;
    
    Sketch::Hash_set minHashes;
    
    const Sketch & sketch = input->sketch;
    int kmerSize = sketch.getKmerSize();
    int mins = sketch.getMinHashesPerWindow();
    int length = input->length;
    char * seq = input->seq;
    float threshold = input->threshold;
    int windowSize = sketch.getWindowSize();
    int best = input->best;
    int selfIndexRef = sketch.getReferenceIndex(input->seqId);
    bool selfMatches = input->selfMatches;
    
    output->seqId = input->seqId;
    
    //cout << "Mins: " << mins << "\t length: " << length << "\tComp: " << compressionFactor << endl;
    
    if ( minusStrand )
    {
        char * seqMinus = new char[length];
        
        for ( int i = 0; i < length; i++ )
        {
            char baseMinus = seq[i];
            
            switch ( baseMinus )
            {
                case 'A': baseMinus = 'T'; break;
                case 'C': baseMinus = 'G'; break;
                case 'G': baseMinus = 'C'; break;
                case 'T': baseMinus = 'A'; break;
                default: break;
            }
            
            seqMinus[length - i - 1] = baseMinus;
        }
        
        //cout << seq << endl;
        seq = seqMinus;
        //cout << seq << endl;
    }
    
    //getMinHashes(minHashes, seq, length, 0, kmerSize, mins);
    
    vector<Sketch::PositionHash> positionHashes;
    
    Sketch::Parameters params;
    
    params.kmerSize = kmerSize;
    params.minHashesPerWindow = mins;
    params.windowed = true;
    params.windowSize = windowSize;
    params.use64 = true;
    
    getMinHashPositions(positionHashes, seq, length, params);
    //
    for ( int i = 0; i < positionHashes.size(); i++ )
    {
        minHashes.insert(positionHashes.at(i).hash);
    }
    
    if ( minusStrand )
    {
        delete [] seq;
    }
    
    // get sorted lists of positions, per reference sequence, that have
    // mutual min-hashes with the query
    //
    PositionsBySequence_umap hits;
    //
    for ( Sketch::Hash_set::const_iterator i = minHashes.begin(); i != minHashes.end(); i++ )
    {
        Sketch::hash_t hash = *i;
        
        if ( sketch.hasLociByHash(hash) )
        {
            const vector<Sketch::Locus> loci = sketch.getLociByHash(hash);
            
            for ( int j = 0; j < loci.size(); j++ )
            {
                const Sketch::Locus & locus = loci.at(j);
                
                if ( verbose ) cout << "Match for hash " << hash << "\t" << locus.sequence << "\t" << locus.position << endl;
                
                if ( locus.sequence != selfIndexRef || selfMatches )
                {
                    hits[locus.sequence].insert(locus.position); // set will be created if needed
                }
            }
        }
    }
    
    for ( PositionsBySequence_umap::iterator i = hits.begin(); i != hits.end(); i++ )
    {
    	using std::set;
    	
        // pointer to the position at the beginning of the window; to be updated
        // as the end of the window is incremented
        //
        set<uint32_t>::const_iterator windowStart = i->second.begin();
        
        if ( verbose ) cout << "Clustering in seq " << i->first << endl;
        
        // the number of positions between the window start and end (inclusive)
        //
        int windowCount = 0;
        
        for ( set<uint32_t>::const_iterator j = i->second.begin(); j != i->second.end(); j++ )
        {
            windowCount++;
            
            if ( verbose ) cout << *windowStart << "\t" << *j << endl;
            
            // update window start if it is too far behind
            //
            while ( windowStart != j && *j > length && *windowStart < *j - length + 1 )
            {
                if ( verbose ) cout << "moving " << *j - length + 1 << endl;
                windowStart++;
                windowCount--;
            }
            
            // extend the right of the window if possible
            //
            while ( j != i->second.end() && *j - *windowStart < length )
            {
                windowCount++;
                j++;
            }
            //
            windowCount--;
            j--;
            
            if ( verbose ) cout << *windowStart << "\t" << *j << endl;
            float score = float(windowCount) / minHashes.size();
            
            if
            (
                score >= threshold &&
                (
                    best == 0 ||
                    output->hits.size() < best ||
                    CommandFind::FindOutput::Hit(i->first, *windowStart, *j, minusStrand, score) < output->hits.top()
                )
            )
            {
                if ( verbose ) cout << input->seqId << '\t' << sketch.getReference(i->first).name << '\t' << *windowStart << '\t' << *j << '\t' << float(windowCount) / mins << endl;
                
                output->hits.push(CommandFind::FindOutput::Hit(i->first, *windowStart, *j, minusStrand, score));
                
                if ( best != 0 && output->hits.size() > best )
                {
                    output->hits.pop();
                }
                
                //break;
                
                for ( set<uint32_t>::const_iterator k = windowStart; k != i->second.end() && *k <= *j; k++ )
                {
                    if ( verbose ) cout << "      " << *k << endl;
                }
            }
        }
    }
    
    //cout << "done\n";
}

bool operator<(const CommandFind::FindOutput::Hit & a, const CommandFind::FindOutput::Hit & b)
{
    // inverse since we want to know the lowest score in our priority queue
    
    if ( a.score == b.score )
    {
        if ( a.ref == b.ref )
        {
            if ( a.start == b.start )
            {
                return b.minusStrand;
            }
            
            return a.start < b.start;
        }
        
        return a.ref < b.ref;
    }
    
    return a.score > b.score;
}

} // namespace mash
