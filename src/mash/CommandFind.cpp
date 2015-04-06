#include "CommandFind.h"
#include "Index.h"
#include <zlib.h>
#include "kseq.h"
#include <iostream>
#include <set>
#include <unordered_set>
#include "ThreadPool.h"

using namespace::std;

KSEQ_INIT(gzFile, gzread)

CommandFind::CommandFind()
: Command()
{
    name = "find";
    description = "Compare query sequences to a reference. <reference> can be a fasta file (gzipped or not) or a mash windowed sketch file (.msw). If it is fasta and a sketch file does not exist (or is out of date), one will be written with the current options and used for future runs if possible. <query> can be fasta or fastq, gzipped or not. Multiple query files can be provided, or \"-\" can be given to read from standard input.";
    argumentString = "<reference> <query> [<query>] ...";
    
    useOption("help");
    useOption("kmer");
    useOption("window");
    useOption("minsWindowed");
    addOption("threshold", Option(Option::Number, "t", "Threshold. This fraction of the query sequence's min-hashes must appear in a query-sized window of a reference sequence for the match to be reported.", "0.2", 0.0, 1.0));
    addOption("threads", Option(Option::Integer, "p", "Parallelism. This many threads will be spawned to perform the find, each one handling on query sequence at a time.", "1"));
    addOption("best", Option(Option::Integer, "b", "Best hit count. This many of the best hits will be reported (0 to report all hits). Score ties are broken by keeping the hit to the earlier reference or to the left-most position.", "0"));
    addOption("self", Option(Option::Boolean, "self", "Allow self matches if query ID appears in reference index.", ""));
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
    bool selfMatches = options.at("self").active;
    
    Index index;
    const string & fileReference = arguments[0];
    
    if ( hasSuffix(fileReference, suffixSketchWindowed) )
    {
    	index.initFromCapnp(fileReference.c_str());
    }
    else
    {
		int kmerSize = options.at("kmer").getArgumentAsNumber();
		int mins = options.at("minsWindowed").getArgumentAsNumber();
		int windowSize = options.at("window").getArgumentAsNumber();
		
		bool indexFileExists = index.initHeaderFromBaseIfValid(fileReference, true);
		
		if
		(
			(options.at("kmer").active && kmerSize != index.getKmerSize()) ||
			(options.at("minsWindowed").active && mins != index.getMinHashesPerWindow()) ||
			(options.at("window").active && windowSize != index.getWindowSize())
		)
		{
			indexFileExists = false;
		}
		
		if ( indexFileExists )
		{
			index.initFromBase(arguments[0], true);
			kmerSize = index.getKmerSize();
			mins = index.getMinHashesPerWindow();
			windowSize = index.getWindowSize();
		}
		else
		{
			vector<string> refArgVector;
			refArgVector.push_back(fileReference);
		
			cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;
			index.initFromSequence(refArgVector, kmerSize, mins, true, windowSize, false);
		
			if ( index.writeToFile() )
			{
				cerr << "Sketch saved for subsequent runs." << endl;
			}
			else
			{
				cerr << "The sketch for " << fileReference << " could not be saved; it will be sketched again next time." << endl;
			}
		}
    }
    
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
        }
        
        gzFile fp = gzdopen(fileno(inStream), "r");
        kseq_t *seq = kseq_init(fp);
        
        while ((l = kseq_read(seq)) >= 0)
        {
            if ( l < index.getKmerSize() )
            {
                continue;
            }
            
            //printf("Query name: %s\tlength: %d\n\n", seq->name.s, l);
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            
            threadPool.runWhenThreadAvailable(new FindInput(index, seq->name.s, seq->seq.s, l, threshold, best, selfMatches));
            
            while ( threadPool.outputAvailable() )
            {
                //cout << "popping\n";
                writeOutput(index, threadPool.popOutputWhenAvailable());
            }
        }
        
        if ( l != -1 )
        {
            printf("ERROR: return value: %d\n", l);
            return 1;
        }
        
        kseq_destroy(seq);
        gzclose(fp);
    }
    
    while ( threadPool.running() )
    {
        writeOutput(index, threadPool.popOutputWhenAvailable());
    }
    
    return 0;
}

void CommandFind::writeOutput(const Index & index, FindOutput * output) const
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
            index.getReference(hit.ref).name << '\t' <<
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
    typedef unordered_map < uint32_t, set<uint32_t> > PositionsBySequence_umap;
    
    bool verbose = false;
    
    Index::Hash_set minHashes;
    
    const Index & index = input->index;
    int kmerSize = index.getKmerSize();
    int mins = index.getMinHashesPerWindow();
    int length = input->length;
    char * seq = input->seq;
    float threshold = input->threshold;
    int windowSize = index.getWindowSize();
    int best = input->best;
    int selfIndexRef = index.getReferenceIndex(input->seqId);
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
    
    vector<Index::PositionHash> positionHashes;
    getMinHashPositions(positionHashes, seq, length, kmerSize, mins, windowSize, 0);
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
    for ( Index::Hash_set::const_iterator i = minHashes.begin(); i != minHashes.end(); i++ )
    {
        Index::hash_t hash = *i;
        //cout << "Hash " << hash << endl;
        
        if ( index.hasLociByHash(hash) )
        {
            const vector<Index::Locus> loci = index.getLociByHash(hash);
            
            for ( int j = 0; j < loci.size(); j++ )
            {
                const Index::Locus & locus = loci.at(j);
                
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
                if ( verbose ) cout << input->seqId << '\t' << index.getReference(i->first).name << '\t' << *windowStart << '\t' << *j << '\t' << float(windowCount) / mins << endl;
                
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
