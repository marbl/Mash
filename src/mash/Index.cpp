#include "Index.h"
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fcntl.h>
#include <map>
#include "kseq.h"
#include "MurmurHash3.h"
#include <assert.h>
#include <queue>
#include <deque>
#include <set>
#include "Command.h" // TEMP for column printing
#include <sys/stat.h>

#define SET_BINARY_MODE(file)
#define CHUNK 16384
KSEQ_INIT(gzFile, gzread)

using namespace std;

typedef map < Index::hash_t, vector<Index::PositionHash> > LociByHash_map;

const vector<Index::Locus> & Index::getLociByHash(hash_t hash) const
{
    return lociByHash.at(hash);
}

int Index::getReferenceIndex(string id) const
{
    if ( referenceIndecesById.count(id) == 1 )
    {
        return referenceIndecesById.at(id);
    }
    else
    {
        return -1;
    }
}

bool Index::initFromBase(const std::string & fileSeq, bool windowed)
{
    file = fileSeq + (windowed ? suffixWindowed : suffix);
    
    struct stat fileInfoSeq;
    struct stat fileInfoIndex;
    
    if ( stat(fileSeq.c_str(), &fileInfoSeq) == -1 )
    {
        return false;
    }
    
    if ( stat(file.c_str(), &fileInfoIndex) == -1 )
    {
        return false;
    }
    
    if ( fileInfoSeq.st_mtime > fileInfoIndex.st_mtime )
    {
        return false;
    }
    
    initFromCapnp(file.c_str());
    return true;
}

int Index::initFromCapnp(const char * file)
{
    // use a pipe to decompress input to Cap'n Proto
    
    int fds[2];
    int piped = pipe(fds);
    
    if ( piped < 0 )
    {
        cerr << "ERROR: could not open pipe for decompression\n";
        return 1;
    }
    
    int forked = fork();
    
    if ( forked < 0 )
    {
        cerr << "ERROR: could not fork for decompression" << endl;
        return 1;
    }
    
    if ( forked == 0 )
    {
        // read from zipped fd and write to pipe
        
        close(fds[0]); // other process's end of pipe
        
        int fd = open(file, O_RDONLY);
        
        if ( fd < 0 )
        {
            cerr << "ERROR: could not open " << file << " for reading." << endl;
            exit(1);
        }
        
        char buffer[1024];
        
        read(fd, buffer, capnpHeaderLength);
        buffer[capnpHeaderLength] = 0;
        
        if ( strcmp(buffer, capnpHeader) != 0 )
        {
            cerr << "ERROR: '" << file << "' does not look like a mash index" << endl;
            exit(1);
        }
        
        int ret = inf(fd, fds[1]);
        if (ret != Z_OK) zerr(ret);
        close(fd);
        exit(ret);
        
        gzFile fileIn = gzopen(file, "rb");
        
        int bytesRead;
        
        // eat header
        //
        gzread(fileIn, buffer, capnpHeaderLength);
        
        printf("header: %s\n", buffer);
        while ( (bytesRead = gzread(fileIn, buffer, sizeof(buffer))) > 0)
        {
            printf("uncompressed: %s\n", buffer);
            write(fds[1], buffer, bytesRead);
        }
        
        gzclose(fileIn);
        close(fds[1]);
        exit(0);
    }
    
    // read from pipe
    
    close(fds[1]); // other process's end of pipe
    
    capnp::ReaderOptions readerOptions;
    
    readerOptions.traversalLimitInWords = 1000000000000;
    readerOptions.nestingLimit = 1000000;
    
    capnp::StreamFdMessageReader message(fds[0], readerOptions);
    capnp::MinHash::Reader reader = message.getRoot<capnp::MinHash>();
    
    capnp::MinHash::ReferenceList::Reader referenceListReader = reader.getReferenceList();
    
    capnp::List<capnp::MinHash::ReferenceList::Reference>::Reader referencesReader = referenceListReader.getReferences();
    references.resize(referencesReader.size());
    
    for ( int i = 0; i < references.size(); i++ )
    {
        capnp::MinHash::ReferenceList::Reference::Reader referenceReader = referencesReader[i];
        
        references[i].name = referenceReader.getName();
        references[i].comment = referenceReader.getComment();
        references[i].length = referenceReader.getLength();
        
        capnp::List<uint64_t>::Reader hashesReader = referenceReader.getHashes();
        
        for ( int j = 0; j < hashesReader.size(); j++ )
        {
            references[i].hashes.insert(hashesReader[j]);
        }
    }
    
    capnp::MinHash::LocusList::Reader locusListReader = reader.getLocusList();
    capnp::List<capnp::MinHash::LocusList::Locus>::Reader lociReader = locusListReader.getLoci();
    
    positionHashesByReference.resize(references.size());
    
    for ( int i = 0; i < lociReader.size(); i++ )
    {
        capnp::MinHash::LocusList::Locus::Reader locusReader = lociReader[i];
        //cout << locusReader.getHash() << '\t' << locusReader.getSequence() << '\t' << locusReader.getPosition() << endl;
        positionHashesByReference[locusReader.getSequence()].push_back(PositionHash(locusReader.getPosition(), locusReader.getHash()));
    }
    
    kmerSize = reader.getKmerSize();
    minHashesPerWindow = reader.getMinHashesPerWindow();
    windowSize = reader.getWindowSize();
    
    //cout << endl << "References:" << endl << endl;
    
    vector< vector<string> > columns(3);
    
    columns[0].push_back("ID");
    columns[1].push_back("Length");
    columns[2].push_back("Name/Comment");
    
    for ( int i = 0; i < references.size(); i++ )
    {
        columns[0].push_back(to_string(i));
        columns[1].push_back(to_string(references[i].length));
        columns[2].push_back(references[i].name + " " + references[i].comment);
    }
    
    //printColumns(columns);
    //cout << endl;
    
    /*
    printf("\nCombined hash table:\n");
    
    cout << "   kmer:  " << kmerSize << endl;
    cout << "   comp:  " << compressionFactor << endl << endl;
    
    for ( LociByHash_umap::iterator i = lociByHash.begin(); i != lociByHash.end(); i++ )
    {
        printf("Hash %u:\n", i->first);
        
        for ( int j = 0; j < i->second.size(); j++ )
        {
            printf("   Seq: %d\tPos: %d\n", i->second.at(j).sequence, i->second.at(j).position);
        }
    }
    
    cout << endl;
    */
    close(fds[0]);
    
    createTables();
    
    return 0;
}

int Index::initFromSequence(const vector<string> & files, int kmerSizeNew, int minHashesPerWindowNew, bool windowedNew, int windowSizeNew, bool concat, int verbosity)
{
    kmerSize = kmerSizeNew;
    minHashesPerWindow = minHashesPerWindowNew;
    windowSize = windowSizeNew;
    windowed = windowedNew;
    
    priority_queue<Index::hash_t> minHashesQueue; // only used for non-windowed
    
    int l;
    int count = 0;
    
    for ( int i = 0; i < files.size(); i++ )
    {
        gzFile fp = gzopen(files[i].c_str(), "r");
        kseq_t *seq = kseq_init(fp);
        
        if ( concat )
        {
            references.resize(references.size() + 1);
            references[references.size() - 1].name = files[i];
        }
        
        while ((l = kseq_read(seq)) >= 0)
        {
            if ( l < kmerSize )
            {
                continue;
            }
            
            if ( windowed )
            {
                positionHashesByReference.resize(count + 1);
            }
            
            if ( ! concat )
            {
                references.resize(references.size() + 1);
            }
            
            if ( verbosity > 0 && windowed ) cout << '>' << seq->name.s << " (" << l << "nt)" << endl << endl;
            //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //printf("seq: %s\n", seq->seq.s);
            //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
            
            Reference & reference = references[references.size() - 1];
            
            if ( ! concat )
            {
                reference.name = seq->name.s;
                
                if ( seq->comment.l > 0 )
                {
                    reference.comment = seq->comment.s;
                }
                
                reference.length = l;
            }
            
            if ( windowed )
            {
                getMinHashPositions(positionHashesByReference[count], seq->seq.s, l, kmerSize, minHashesPerWindow, windowSize, verbosity);
            }
            else
            {
                addMinHashes(reference.hashes, minHashesQueue, seq->seq.s, l, kmerSize, minHashesPerWindow);
            }
            
            if ( ! concat )
            {
                count++;
            }
        }
        
        if ( l != -1 )
        {
            printf("ERROR: return value: %d\n", l);
            return 1;
        }
        
        if ( concat )
        {
            count++;
        }
        
        kseq_destroy(seq);
        gzclose(fp);
    }
    /*
    printf("\nCombined hash table:\n\n");
    
    for ( LociByHash_umap::iterator i = lociByHash.begin(); i != lociByHash.end(); i++ )
    {
        printf("Hash %u:\n", i->first);
        
        for ( int j = 0; j < i->second.size(); j++ )
        {
            printf("   Seq: %d\tPos: %d\n", i->second.at(j).sequence, i->second.at(j).position);
        }
    }
    */
    
    createTables();
    
    return 0;
}

bool Index::writeToFile() const
{
    return writeToCapnp(file.c_str()) == 0;
}

int Index::writeToCapnp(const char * file) const
{
    // use a pipe to compress Cap'n Proto output
    
    int fds[2];
    int piped = pipe(fds);
    
    if ( piped < 0 )
    {
        cerr << "ERROR: could not open pipe for compression\n";
        return 1;
    }
    
    int forked = fork();
    
    if ( forked < 0 )
    {
        cerr << "ERROR: could not fork for compression\n";
        return 1;
    }
    
    if ( forked == 0 )
    {
        // read from pipe and write to compressed file
        
        close(fds[1]); // other process's end of pipe
        
        int fd = open(file, O_CREAT | O_WRONLY | O_TRUNC, 0644);
        
        if ( fd < 0 )
        {
            cerr << "ERROR: could not open " << file << " for writing.\n";
            exit(1);
        }
        
        // write header
        //
        write(fd, capnpHeader, capnpHeaderLength);
        
        int ret = def(fds[0], fd, Z_DEFAULT_COMPRESSION);
        if (ret != Z_OK) zerr(ret);
        exit(ret);
        
        char buffer[1024];
        gzFile fileOut = gzopen(file, "ab");
        
        int bytesRead;
        
        while ( (bytesRead = read(fds[0], buffer, sizeof(buffer))) > 0)
        {
            printf("compressing: %s\n", buffer);
            gzwrite(fileOut, buffer, bytesRead);
        }
        
        gzclose(fileOut);
        close(fds[0]);
        exit(0);
    }
    
    // write to pipe
    
    close(fds[0]); // other process's end of pipe
    
    capnp::MallocMessageBuilder message;
    capnp::MinHash::Builder builder = message.initRoot<capnp::MinHash>();
    
    capnp::MinHash::ReferenceList::Builder referenceListBuilder = builder.initReferenceList();
    
    capnp::List<capnp::MinHash::ReferenceList::Reference>::Builder referencesBuilder = referenceListBuilder.initReferences(references.size());
    
    for ( int i = 0; i < references.size(); i++ )
    {
        capnp::MinHash::ReferenceList::Reference::Builder referenceBuilder = referencesBuilder[i];
        
        referenceBuilder.setName(references[i].name);
        referenceBuilder.setComment(references[i].comment);
        referenceBuilder.setLength(references[i].length);
        
        if ( references[i].hashes.size() != 0 )
        {
            const Hash_set & hashes = references[i].hashes;
            capnp::List<uint64_t>::Builder hashesBuilder = referenceBuilder.initHashes(hashes.size());
            
            int index = 0;
            
            for ( Hash_set::const_iterator j = hashes.begin(); j != hashes.end(); j++ )
            {
                hashesBuilder.set(index, *j);
                index++;
            }
        }
    }
    
    int locusCount = 0;
    
    for ( int i = 0; i < positionHashesByReference.size(); i++ )
    {
        locusCount += positionHashesByReference.at(i).size();
    }
    
    capnp::MinHash::LocusList::Builder locusListBuilder = builder.initLocusList();
    capnp::List<capnp::MinHash::LocusList::Locus>::Builder lociBuilder = locusListBuilder.initLoci(locusCount);
    
    int locusIndex = 0;
    
    for ( int i = 0; i < positionHashesByReference.size(); i++ )
    {
        for ( int j = 0; j < positionHashesByReference.at(i).size(); j++ )
        {
            capnp::MinHash::LocusList::Locus::Builder locusBuilder = lociBuilder[locusIndex];
            locusIndex++;
            
            locusBuilder.setSequence(i);
            locusBuilder.setPosition(positionHashesByReference.at(i).at(j).position);
            locusBuilder.setHash(positionHashesByReference.at(i).at(j).hash);
        }
    }
    
    builder.setKmerSize(kmerSize);
    builder.setMinHashesPerWindow(minHashesPerWindow);
    builder.setWindowSize(windowSize);
    
    writeMessageToFd(fds[1], message);
    close(fds[1]);
    
    return 0; // TODO
}

void Index::createTables()
{
    for ( int i = 0; i < references.size(); i++ )
    {
        referenceIndecesById[references[i].name] = i;
    }
    
    for ( int i = 0; i < positionHashesByReference.size(); i++ )
    {
        for ( int j = 0; j < positionHashesByReference.at(i).size(); j++ )
        {
        	const PositionHash & positionHash = positionHashesByReference.at(i).at(j);
        	
	        lociByHash[positionHash.hash].push_back(Locus(i, positionHash.position));
        }
    }
}

void addMinHashes(Index::Hash_set & minHashes, priority_queue<Index::hash_t> & minHashesQueue, char * seq, uint32_t length, int kmerSize, int mins)
{
    //cout << "mins: " << mins << endl << endl;
    
    // uppercase
    //
    for ( int i = 0; i < length; i++ )
    {
        if ( seq[i] > 90 )
        {
            seq[i] -= 32;
        }
    }
    
    for ( int i = 0; i < length - kmerSize + 1; i++ )
    {
        // repeatedly skip kmers with bad characters
        //
        for ( int j = i; j < i + kmerSize && i + kmerSize <= length; j++ )
        {
            char c = seq[j];
            
            if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' )
            {
                i = j + 1; // skip to past the bad character
                break;
            }
        }
        
        if ( i + kmerSize > length )
        {
            // skipped to end
            break;
        }
        
        Index::hash_t hash = getHash(seq + i, kmerSize);
        
        //if ( i % 1000000 == 0 )
        {
            //printf("   At position %d\n", i);
        }
        
        if
        (
            (
                minHashesQueue.size() < mins ||
                hash < minHashesQueue.top()
            )
            && minHashes.count(hash) == 0
        )
        {
            minHashes.insert(hash);
            minHashesQueue.push(hash);
            
            if ( minHashesQueue.size() > mins )
            {
                minHashes.erase(minHashesQueue.top());
                minHashesQueue.pop();
            }
        }
    }
}

Index::hash_t getHash(const char * seq, int length)
{
    Index::hash_t hash = 0;
    
#ifdef ARCH_32
    uint32_t data[4];
    MurmurHash3_x86_32(seq, length, seed, data);
    hash = data[0];
#else
    uint64_t data[2];
    MurmurHash3_x64_128(seq, length, seed, data);
    hash = data[0];
#endif    

    return hash;
}

void getMinHashPositions(vector<Index::PositionHash> & positionHashes, char * seq, uint32_t length, int kmerSize, int minHashesPerWindow, int windowSize, int verbosity)
{
    int mins = minHashesPerWindow;
    int nextValidKmer = 0;
    
    struct CandidateLocus
    {
        CandidateLocus(int positionNew)
            :
            position(positionNew),
            isMinmer(false)
            {}
        
        int position;
        bool isMinmer;
    };
    
    if ( windowSize > length - kmerSize + 1 )
    {
        windowSize = length - kmerSize + 1;
    }
    
    if ( verbosity > 1 ) cout << seq << endl << endl;
    map<Index::hash_t, deque<CandidateLocus>> candidatesByHash;
    
    queue<map<Index::hash_t, deque<CandidateLocus>>::iterator> windowQueue;
    map<Index::hash_t, deque<CandidateLocus>>::iterator maxMinmer = candidatesByHash.end();
    map<Index::hash_t, deque<CandidateLocus>>::iterator newCandidates;
    
    int unique = 0;
    
    for ( int i = 0; i < length - kmerSize + 1; i++ )
    {
        if ( i >= nextValidKmer )
        {
            for ( int j = i; j < i + kmerSize; j++ )
            {
                char c = seq[j];
                
                if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' )
                {
                    //nextValidKmer = j + 1;
                    break;
                }
            }
        }
        
        if ( i < nextValidKmer && verbosity > 1 )
        {
            cout << "  [";
        
            for ( int j = i; j < i + kmerSize; j++ )
            {
                cout << seq[j];
            }
            
            cout << "]" << endl;
        }
        
        if ( i >= nextValidKmer )
        {
            Index::hash_t hash = getHash(seq + i, kmerSize);
            
            if ( verbosity > 1 )
            {
                cout << "   ";
            
                for ( int j = i; j < i + kmerSize; j++ )
                {
                    cout << seq[j]; 
                }
            
                cout << "   " << i << '\t' << hash << endl;
            }
            
            pair<map<Index::hash_t, deque<CandidateLocus>>::iterator, bool> inserted =
                candidatesByHash.insert(pair<Index::hash_t, deque<CandidateLocus>>(hash, deque<CandidateLocus>()));
            newCandidates = inserted.first;
            newCandidates->second.push_back(CandidateLocus(i));
            
            if
            (
                inserted.second && // inserted; decrement maxMinmer if...
                (
                    (
                        // ...just reached number of mins
                        
                        maxMinmer == candidatesByHash.end() &&
                        candidatesByHash.size() == mins
                    ) ||
                    (
                        // ...inserted before maxMinmer
                        
                        maxMinmer != candidatesByHash.end() &&
                        newCandidates->first < maxMinmer->first
                    )
                )
            )
            {
                maxMinmer--;
                
                if ( i >= windowSize )
                {
                    unique++;
                }
            }
        }
        else
        {
            newCandidates = candidatesByHash.end();
        }
        
        windowQueue.push(newCandidates);
        map<Index::hash_t, deque<CandidateLocus>>::iterator windowFront = candidatesByHash.end();
        
        if ( windowQueue.size() > windowSize )
        {
            windowFront = windowQueue.front();
            windowQueue.pop();
            
            if ( verbosity > 1 ) cout << "   \tPOP: " << windowFront->first << endl;
        }
        
        if ( windowFront != candidatesByHash.end() )
        {
            deque<CandidateLocus> & frontCandidates = windowFront->second;
            
            if ( frontCandidates.front().isMinmer )
            {
                if ( verbosity > 1 ) cout << "   \t   minmer: " << frontCandidates.front().position << '\t' << windowFront->first << endl;
                positionHashes.push_back(Index::PositionHash(frontCandidates.front().position, windowFront->first));
            }
            
            if ( frontCandidates.size() > 1 )
            {
                frontCandidates.pop_front();
            
                if ( maxMinmer == candidatesByHash.end() || ( i >= windowSize && windowFront->first <= maxMinmer->first) )
                {
                    frontCandidates.front().isMinmer = true;
                }
            }
            else
            {
                if ( maxMinmer != candidatesByHash.end() && windowFront->first <= maxMinmer->first )
                {
                    maxMinmer++;
                    
                    if ( maxMinmer != candidatesByHash.end() )
                    {
                        maxMinmer->second.front().isMinmer = true;
                    }
                    
                    unique++;
                }
            
                candidatesByHash.erase(windowFront);
            }
        }
        
        if ( i == windowSize - 1 )
        {
            // first complete window; mark minmers
            
            for ( map<Index::hash_t, deque<CandidateLocus>>::iterator j = candidatesByHash.begin(); j != maxMinmer; j++ )
            {
                j->second.front().isMinmer = true;
            }
            
            if ( maxMinmer != candidatesByHash.end() )
            {
                maxMinmer->second.front().isMinmer = true;
            }
            
            unique++;
        }
        
        if ( newCandidates != candidatesByHash.end() && i >= windowSize && (maxMinmer == candidatesByHash.end() || newCandidates->first <= maxMinmer->first) )
        {
            newCandidates->second.front().isMinmer = true;
        }
        
        if ( verbosity > 1 )
        {
            for ( map<Index::hash_t, deque<CandidateLocus>>::iterator j = candidatesByHash.begin(); j != candidatesByHash.end(); j++ )
            {
                cout << "   \t" << j->first;
                
                if ( j == maxMinmer )
                {
                     cout << "*";
                }
                
                for ( deque<CandidateLocus>::iterator k = j->second.begin(); k != j->second.end(); k++ )
                {
                    cout << '\t' << k->position;
                    
                    if ( k->isMinmer )
                    {
                        cout << '!';
                    }
                }
                
                cout << endl;
            }
        }
    }
    
    // finalize remaining min-hashes from the last window
    //
    while ( windowQueue.size() > 0 )
    {
        map<Index::hash_t, deque<CandidateLocus>>::iterator windowFront = windowQueue.front();
        windowQueue.pop();
        
        if ( windowFront != candidatesByHash.end() )
        {
            deque<CandidateLocus> & frontCandidates = windowFront->second;
            
            if ( frontCandidates.size() > 0 )
            {
                if ( frontCandidates.front().isMinmer )
                {
                    if ( verbosity > 1 ) cout << "   \t   minmer:" << frontCandidates.front().position << '\t' << windowFront->first << endl;
                    positionHashes.push_back(Index::PositionHash(frontCandidates.front().position, windowFront->first));
                }
                
                frontCandidates.pop_front();
            }
        }
    }
    
    if ( verbosity > 1 )
    {
        cout << endl << "Minmers:" << endl;
    
        for ( int i = 0; i < positionHashes.size(); i++ )
        {
            cout << "   " << positionHashes.at(i).position << '\t' << positionHashes.at(i).hash << endl;
        }
        
        cout << endl;
    }
    
    if ( verbosity > 0 ) cout << "   " << positionHashes.size() << " minmers across " << length - windowSize - kmerSize + 2 << " windows (" << unique << " windows with distinct minmer sets)." << endl << endl;
}

// The following functions are adapted from http://www.zlib.net/zpipe.c


/* Compress from file source to file dest until EOF on source.
   def() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_STREAM_ERROR if an invalid compression
   level is supplied, Z_VERSION_ERROR if the version of zlib.h and the
   version of the library linked do not match, or Z_ERRNO if there is
   an error reading or writing the files. */
int def(int fdSource, int fdDest, int level)
{
    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        return ret;

    /* compress until end of file */
    do {
        strm.avail_in = read(fdSource, in, CHUNK);
        if (strm.avail_in == -1) {
            (void)deflateEnd(&strm);
            return Z_ERRNO;
        }
        flush = strm.avail_in == 0 ? Z_FINISH : Z_NO_FLUSH;
        strm.next_in = in;

        /* run deflate() on input until output buffer not full, finish
           compression if all of source has been read in */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = deflate(&strm, flush);    /* no bad return value */
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            have = CHUNK - strm.avail_out;
            if (write(fdDest, out, have) != have) {
                (void)deflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);
        assert(strm.avail_in == 0);     /* all input will be used */

        /* done when last data in file processed */
    } while (flush != Z_FINISH);
    assert(ret == Z_STREAM_END);        /* stream will be complete */

    /* clean up and return */
    (void)deflateEnd(&strm);
    return Z_OK;
}

/* Decompress from file source to file dest until stream ends or EOF.
   inf() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_DATA_ERROR if the deflate data is
   invalid or incomplete, Z_VERSION_ERROR if the version of zlib.h and
   the version of the library linked do not match, or Z_ERRNO if there
   is an error reading or writing the files. */
int inf(int fdSource, int fdDest)
{
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return ret;

    /* decompress until deflate stream ends or end of file */
    do {
        strm.avail_in = read(fdSource, in, CHUNK);
        if (strm.avail_in == -1) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            if (write(fdDest, out, have) != have) {
                (void)inflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

/* report a zlib or i/o error */
void zerr(int ret)
{
    fputs("zpipe: ", stderr);
    switch (ret) {
    case Z_ERRNO:
        if (ferror(stdin))
            fputs("error reading stdin\n", stderr);
        if (ferror(stdout))
            fputs("error writing stdout\n", stderr);
        break;
    case Z_STREAM_ERROR:
        fputs("invalid compression level\n", stderr);
        break;
    case Z_DATA_ERROR:
        fputs("invalid or incomplete deflate data\n", stderr);
        break;
    case Z_MEM_ERROR:
        fputs("out of memory\n", stderr);
        break;
    case Z_VERSION_ERROR:
        fputs("zlib version mismatch!\n", stderr);
    }
}
