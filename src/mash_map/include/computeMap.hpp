/**
 * @file    computeMap.hpp
 * @brief   implments the sequence mapping logic
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_MAP_HPP 
#define SKETCH_MAP_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <zlib.h>  

//Own includes
#include "map_parameters.hpp"
#include "commonFunc.hpp"
#include "winSketch.hpp"
#include "map_stats.hpp"

//External includes


namespace skch
{
  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages
   */
  class Map
  {
    //private members
    const std::vector<std::string> &queryFiles;
    const std::string &outputFile;

    int kmerSize;
    int windowSize;
    int minimumReadLength;

    //Threshold for identity
    float percentageIdentity;

    //alphabet size
    int alphabetSize;

    //reference sketch
    const skch::Sketch &refSketch;

    //(seq. id, start ref. position, end position) for L1
    typedef std::tuple<seqno_t, offset_t, offset_t> candidateLocus_t;

    //Container type for saving read sketches during L1 and L2 both
    typedef Sketch::MI_Type MinVec_Type;

    typedef Sketch::MIIter_t MIIter_t;

    //(seq. id, begin iterator over minimizerIndex, end iterator over minimizerIndex, 
    //  equal minhash count) for L2
    //  [begin iterator, end iterator) represents the mapped region on reference
    typedef std::tuple<seqno_t, MIIter_t, MIIter_t, int> mapLocus_t;

    //Make the default constructor private, non-accessible
    Map();

    const int NA = -1;    //hash 'Not Available' marker

    int stepDistance;

    public:

    /**
     * @brief   constructor
     */
    Map(const skch::Parameters &param, const skch::Sketch &refsketch) :
      queryFiles(param.querySequences),
      outputFile(param.outFileName),
      percentageIdentity(param.percentageIdentity), 
      minimumReadLength(param.minReadLength),
      alphabetSize(param.alphabetSize),
      refSketch(refsketch) 
    {
      //Sketching parameters should match the reference for correct jaccard estimate
      kmerSize = refsketch.getKmerSize();
      windowSize = refsketch.getWindowSize();

      //Keeping step sizes to be 10 times window size
      //Example, 1000 step size for win size 100
      this->stepDistance = windowSize * 10;

      this->mapQuery();
    }

    private:

    /**
     * @brief   parse over sequences in query file and map each on the reference
     */
    void mapQuery()
    {
      //Count of reads mapped by us
      //Some reads are dropped because of short length
      seqno_t totalReadsPickedForMapping = 0;
      seqno_t totalReadsMapped = 0;
      seqno_t seqCounter = 0;

      std::ofstream outstrm(this->outputFile);

      for(const auto &fileName : this->queryFiles)
      {
        //Open the file using kseq
        FILE *file = fopen(fileName.c_str(), "r");
        gzFile fp = gzdopen(fileno(file), "r");
        kseq_t *seq = kseq_init(fp);


#ifdef DEBUG
        std::cout << "INFO, skch::Map::mapQuery, mapping reads in " << fileName << std::endl;
#endif

        //size of sequence
        offset_t len;

        while ((len = kseq_read(seq)) >= 0) 
        {
          //Is the read too short?
          if(len < this->windowSize || len < this->kmerSize || len < this->minimumReadLength)
          {
            seqCounter++;

#ifdef DEBUG
            std::cout << "WARNING, skch::Map::mapQuery, read is not long enough for mapping" << std::endl;
#endif

            continue;
          }
          else 
          {
            //Minimizers in the query
            MinVec_Type minimizerTableQuery;

            //L1 Mapping
            std::vector<candidateLocus_t> l1Mappings; 
            doL1Mapping(seq, seqCounter, minimizerTableQuery, l1Mappings);

            //L2 Mapping
            if ( doL2Mapping(seq, seqCounter, outstrm, minimizerTableQuery, l1Mappings) )
              totalReadsMapped++;

            totalReadsPickedForMapping++;
            seqCounter++;
          }
        }

        //Close the input file
        kseq_destroy(seq);  
        gzclose(fp);  
      }

      outstrm.close();

      std::cout << "INFO, skch::Map::mapQuery, [count of mapped reads, reads qualified for mapping, total input reads] = [" << totalReadsMapped << ", " << totalReadsPickedForMapping << ", " << seqCounter << "]" << std::endl;

    }

    /**
     * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
     * @details     The count of hits that should occur within a region on the reference is 
     *              determined by the threshold similarity
     *              The resulting start and end target offsets on reference is (are) an 
     *              overestimate of the mapped region. Computing better bounds is left for
     *              the following L2 stage.
     * @param[in]   seq         kseq fasta/q parser
     * @param[out]  l1Mappings  all the read mapping locations
     */
    template <typename KSEQ, typename Vec>
    void doL1Mapping(KSEQ *seq, seqno_t seqCounter, MinVec_Type &minimizerTableQuery, Vec &l1Mappings)
    {
      //Vector of positions of all the hits 
      std::vector< std::pair<seqno_t, offset_t> > seedHitsL1;

      ///1. Compute the minimizers

      CommonFunc::addMinimizers(minimizerTableQuery, seq, this->kmerSize, this->windowSize, alphabetSize);

#ifdef DEBUG
      std::cout << "INFO, skch::Map:doL1Mapping, read id " << seqCounter << ", minimizer count = " << minimizerTableQuery.size() << "\n";
#endif

      ///2. Find the hits in the reference, pick 's' unique minimizers as seeds, 
      
      std::sort(minimizerTableQuery.begin(), minimizerTableQuery.end(), CommonFunc::TpleComp<0>());

      //note : unique preserves the original relative order of elements 
      auto uniqEndIter = std::unique(minimizerTableQuery.begin(), minimizerTableQuery.end(), CommonFunc::TpleComp<0, std::equal_to>());

      //This is the sketch size for estimating jaccard
      int s = std::distance(minimizerTableQuery.begin(), uniqEndIter);

      int totalMinimizersPicked = 0;

      for(auto it = minimizerTableQuery.begin(); it != uniqEndIter; it++)
      {
        //Check if hash value exists in the reference lookup index
        auto seedFind = refSketch.minimizerPosLookupIndex.find(std::get<0>(*it));

        if(seedFind != refSketch.minimizerPosLookupIndex.end())
        {
          auto hitPositionList = seedFind->second;

          //Save the positions (Ignore high frequency hits)
          if(hitPositionList.size() < refSketch.getFreqThreshold())
          {
            seedHitsL1.insert(seedHitsL1.end(), hitPositionList.begin(), hitPositionList.end());
          }

        }
      }

#ifdef DEBUG
      std::cout << "INFO, skch::Map:doL1Mapping, read id " << seqCounter << ", Count of L1 hits in the reference = " << seedHitsL1.size() << "\n";
#endif

      int minimumHits = Stat::estimateMinimumHitsRelaxed(s, this->kmerSize, this->percentageIdentity);


      this->computeL1CandidateRegions(seedHitsL1, minimumHits, seq->seq.l, l1Mappings);

#ifdef DEBUG
      std::cout << "INFO, skch::Map:doL1Mapping, read id " << seqCounter << ", minimum hits required for a candidate = " << minimumHits << "\n";
      std::cout << "INFO, skch::Map:doL1Mapping, read id " << seqCounter << ", Count of L1 candidate regions = " << l1Mappings.size() << "\n";
#endif


#ifdef DEBUG
      for(auto &e : l1Mappings)
        std::cout << "INFO, skch::Map:doL1Mapping, read id " << seqCounter << ", L1 candidate : [" << this->refSketch.metadata[std::get<0>(e)].first << " " << this->refSketch.metadata[std::get<0>(e)].second << " " << std::get<1>(e) << " " << std::get<2>(e) << "]\n";
#endif

    }

    /**
     * @brief       Helper function to doL1Mapping()
     * @param[out]  l1Mappings  all the read mapping locations
     */
    template <typename Vec1, typename Vec2>
    void computeL1CandidateRegions(Vec1 &seedHitsL1, int minimumHits, offset_t len, Vec2 &l1Mappings)
    {
      if(minimumHits < 1)
        minimumHits = 1;

      //Sort all the hit positions
      std::sort(seedHitsL1.begin(), seedHitsL1.end());

      for(auto it = seedHitsL1.begin(); it != seedHitsL1.end(); it++)
      {
        if(std::distance(it, seedHitsL1.end()) >= minimumHits)
        {
          auto it2 = it + minimumHits -1;
          //[it .. it2] are 'minimumHits' consecutive hits 
          
          //Check if consecutive hits are close enough
          if(it2->first == it->first && it2->second - it->second < len)
          {
            //Save <1st pos --- 2nd pos>
            candidateLocus_t candidate(it->first, std::max(0, it2->second - len), it->second + len);

            //Check if this candidate overlaps with last inserted one
            auto lst = l1Mappings.end(); lst--;

            //match seq_no and see if this candidate begins before last element ends
            if(l1Mappings.size() > 0 && std::get<0>(candidate) == std::get<0>(*lst) && std::get<2>(*lst) >= std::get<1>(candidate))
            {
              //Push the end pos of last candidate locus further out
              std::get<2>(*lst) = std::get<2>(candidate);
            }
            else
              l1Mappings.push_back(candidate);
          }
        }
      }
    }

    /**
     * @brief                   Revise L1 candidate regions to more precise locations
     * @param[in]   seq         kseq fasta/q parser
     * @param[in]   l1Mappings  candidate regions for query sequence found at L1
     * @return      T/F         True if atleast 1 mapping region is proposed
     */
    template <typename KSEQ, typename Vec>
    bool doL2Mapping(KSEQ *seq, seqno_t seqCounter, std::ofstream &outstrm, MinVec_Type &minimizerTableQuery, 
                    Vec &l1Mappings)
    {
      ///1. Compute minimum s unique hashes among minimizers 
      //Sort by minimizers' hash values
      if(!std::is_sorted(minimizerTableQuery.begin(), minimizerTableQuery.end(), CommonFunc::TpleComp<0>()))
          std::sort(minimizerTableQuery.begin(), minimizerTableQuery.end(), CommonFunc::TpleComp<0>());

      //compute unique minimizers, note that unique preserves the original relative order of unique elements  
      auto uniqEndIter = std::unique(minimizerTableQuery.begin(), minimizerTableQuery.end(), CommonFunc::TpleComp<0, std::equal_to>());

      //This is the sketch size for estimating jaccard
      int s = std::distance(minimizerTableQuery.begin(), uniqEndIter);

      /*
       * Type of container used for walking read over the reference superwindow
       * We choose a std map : 
       *    [hash value -> (location in the reference superwindow, offset in query)]
       * map preserves the sorted order of hashes, so computing jaccard similarity 
       * is a linear walk, each time for the sliding window
       */
      typedef std::map< hash_t, std::pair<offset_t, offset_t> > slidingMapType;

      slidingMapType slidingWindowMinhashes;

      for(auto it = minimizerTableQuery.begin(); it != uniqEndIter; it++)
        slidingWindowMinhashes.emplace(std::get<0>(*it), std::make_pair(NA, std::get<2>(*it)));     //[hash value] -> (NA, offset in query)

      bool mappingReported = false;

      ///2. Walk the read over the candidate regions and compute the jaccard similarity with minimum s sketches
      for(auto &candidateLocus: l1Mappings)
      {
        //Reset the map to just have the query minimizers
        slidingMapType slidingWindowMinhashesCpy = slidingWindowMinhashes;

        mapLocus_t l2;
        computeL2MappedRegions(slidingWindowMinhashesCpy, s, seq->seq.l, candidateLocus, l2);

        //Compute mash distance using calculated jaccard
        float mash_dist = Stat::j2md(1.0 * std::get<3>(l2)/s, this->kmerSize);

        //Compute lower bound to mash distance within 90% confidence interval
        float mash_dist_lower_bound = Stat::md_lower_bound(mash_dist, s, this->kmerSize, 0.9);

        float nucIdentity = 100 * (1 - mash_dist);
        float nucIdentityUpperBound = 100 * (1 - mash_dist_lower_bound);

        /*
         * Compute addtional statistics to filter out false mappings
         * 1. Unique Minimizers in the reference
         */
        std::vector<float> mappingStatistics;

        slidingWindowMinhashesCpy = slidingWindowMinhashes;
        computeMappingStatistics(slidingWindowMinhashesCpy, std::get<1>(l2),std::get<2>(l2), 
            s, seq->seq.l,
            mappingStatistics); 

        /*An alignment is reported if 
         *    the nucleotide identity is >= the percentage identity threshold
         *    estimated complexity (uniqueness of minimizers) is above 75%
         */
        if(nucIdentityUpperBound >= this->percentageIdentity && mappingStatistics[0] >= 0.75)
        {
          outstrm << seq->name.s << " " << seq->seq.l 
            << " 0 " << seq->seq.l - 1 << " +/- " 
            << this->refSketch.metadata[std::get<0>(l2)].first 
            << " " << this->refSketch.metadata[std::get<0>(l2)].second 
            << " " << std::get<2>(*std::get<1>(l2)) << " " << std::get<2>(*std::get<2>(l2))  
            << " " << nucIdentity;

            //Print some statistics
            //minhash matching count, minhash s value
            outstrm << " " << std::get<3>(l2) << " " << s << " " << nucIdentityUpperBound;

            //statistics vector
            for(auto &e : mappingStatistics)
              outstrm << " " << e;
          
          outstrm << "\n";

          mappingReported = true;
        }

      }


      return mappingReported;
    }

    /**
     * @brief       Helper function to doL2Mapping()
     * @param[in]   slidingWindowMinhashes    container for computing min s hashes along the sliding window of size len
     * @param[in]   s                         count of hashes over which jaccard similarity is estimated
     * @param[in]   len                       length of query (read) sequence
     * @param[in]   candidateLocus            candidate region computed at L1 stage
     * @param[out]  l2Mappings                best read mapping location within the candidateLocus
     */
    template <typename MapT>
      void computeL2MappedRegions(MapT &slidingWindowMinhashes, int s, offset_t len, candidateLocus_t &candidateLocus, mapLocus_t &l2Mapping)
      {
        //Sequence # in the reference
        seqno_t refSequenceId = std::get<0>(candidateLocus);

        //Start of the candidate offset in the reference sequence
        offset_t currentStart = std::get<1>(candidateLocus);


        /*
         * [currentStart, currentStart + len - 1] will be the first reference 
         * superwindow against which we compute the jaccard similarity of query
         */

        auto startPosition = std::make_pair(refSequenceId, currentStart);
        auto lastPosition = std::make_pair(refSequenceId, currentStart + len - 1);

#ifdef DEBUG
        if(currentStart + len - 1 > std::get<2>(candidateLocus))
          std::cout << "WARNING, skch::Map::computeL2MappedRegions, candidate region shorter than read length" << std::endl;
#endif

        /*
         * [l_iter, u_iter) represents the superwindow iterator range of minimizers 
         * in the reference sketch index against which jaccard similarity of the 
         * query is estimated
         *
         * After each window shift, we advance the range [l_iter, u_iter) forward
         * depending on the 'stepDistance'
         *
         * Therefore, each shift requires 
         *      1. deletion of [l_iter (previous), l_iter) elements from the window
         *      2. addition of [u_iter (previous), u_iter) elements into the window
         */

        auto iterRange = this->refSketch.getIndexRange(startPosition, lastPosition);
        auto l_iter = iterRange.first;
        auto u_iter = iterRange.second;

        //Save the best mapping locations here
        MIIter_t bestStartIter = l_iter, bestEndIter = u_iter;
        int maxSharedMinimizers = 0;

        //Previous range is null at the beginning
        auto l_prevIter = l_iter;
        auto u_prevIter = l_iter;

        //Keep sliding till query crosses the candidate's end position
        //Check if u_iter offset is less than L1 candidate boundary, and reference sequence id should match
        while(std::get<2>(*u_iter) < std::get<2>(candidateLocus) && 
            std::get<1>(*u_iter) == std::get<0>(candidateLocus))
        {
          assert(std::distance(l_iter , u_iter) > 0);

          /// 1. Push minimizers in the [u_prevIter, u_iter) range to the window

          for(auto it = u_prevIter; it != u_iter; it++)
          {
            hash_t hashvalue = std::get<0>(*it);
            offset_t offsetinReference = std::get<2>(*it);

            //if hash doesn't exist in window
            if(slidingWindowMinhashes.find(hashvalue) == slidingWindowMinhashes.end())
              slidingWindowMinhashes[hashvalue] = std::make_pair(offsetinReference, NA);   //add the hash to window
            else
              slidingWindowMinhashes[hashvalue].first = offsetinReference;                 //just revise the offset 
          }

          /// 2. Remove minimizers in the [l_prevIter, l_iter) range from the window 

          for(auto it = l_prevIter; it != l_iter; it++)
          {
            hash_t hashvalue = std::get<0>(*it);
            offset_t offsetinReference = std::get<2>(*it);

            //hash must exist in the map, by our logic
            assert(slidingWindowMinhashes.find(hashvalue) != slidingWindowMinhashes.end());

            //only if hash and offset match
            if(slidingWindowMinhashes[hashvalue].first == offsetinReference) 
            {
              //if hash belongs to query, then update the reference offset to NA
              if(slidingWindowMinhashes[hashvalue].second != NA)
                slidingWindowMinhashes[hashvalue].first = NA;
              else
              {
                //erase hash from map
                auto delIter = slidingWindowMinhashes.find(hashvalue);
                slidingWindowMinhashes.erase(delIter);
              }
            }
          }

          assert(slidingWindowMinhashes.size() >= s);

          /// 3. Compute jaccard similarity now

          int uniqueHashes = 0, currentSharedMinimizers = 0;

          for (auto it = slidingWindowMinhashes.cbegin(); it != slidingWindowMinhashes.cend(); ++it)
          {
            auto offsetValues = it->second;

            if(offsetValues.first != NA && offsetValues.second != NA)
              currentSharedMinimizers += 1;

            uniqueHashes++;
            if(uniqueHashes == s)
              break;
          }

          //Is it the best superwindow so far?
          if(currentSharedMinimizers > maxSharedMinimizers)
          {
            bestStartIter = l_iter;
            bestEndIter = u_iter;
            maxSharedMinimizers = currentSharedMinimizers;
          }

          /// 4. Advance the iterator range

          //Shift the window to right by stepDistance

          //Backup the previous iterator positions 
          l_prevIter = l_iter; u_prevIter = u_iter;

          //Compute the new iterator positions
          if (!this->refSketch.advanceIter(l_iter, 20) || !this->refSketch.advanceIter(u_iter, 20))
            break;
        }

        //Save the result
        l2Mapping = std::make_tuple(refSequenceId, bestStartIter, bestEndIter, maxSharedMinimizers);
      }

    /**
     * @brief                               compute additional statistics for mapping (post-L2)
     * @param[in]   slidingWindowMinhashes    
     * @param[in]   refSequenceId           begin iterator on reference index by L2 mapping
     * @param[in]   refEndItr               end iterator on reference index by L2 mapping 
     * @param[in]   len                     length of query sequence
     * @param[in]   s
     * @param[out]  statValues     
     */
    template <typename MapT, typename statVec_t>
      void computeMappingStatistics(MapT &slidingWindowMinhashes, MIIter_t refStartItr, MIIter_t refEndItr,
          int s, offset_t len,
          statVec_t &statValues)
      {
        //Push minimizers in the [refStartItr, refEndItr) range to the map
        for(auto it = refStartItr; it != refEndItr; it++)
        {
          hash_t hashvalue = std::get<0>(*it);
          offset_t offsetinReference = std::get<2>(*it);

          //if hash doesn't exist in window
          if(slidingWindowMinhashes.find(hashvalue) == slidingWindowMinhashes.end())
            slidingWindowMinhashes[hashvalue] = std::make_pair(offsetinReference, NA);   //add the hash to window
          else
            slidingWindowMinhashes[hashvalue].first = offsetinReference;                 //just revise the offset 
        }

        ///1. Unique Minimizers in the referece (complexity)
        {
          //Indicator of complexity in the mapped region
          //range (0 - 1.x], lower value means lower complexity in the region
          float referenceDNAComplexity; 

          uint64_t uniqueHashCount = 0;

          for (auto it = slidingWindowMinhashes.cbegin(); it != slidingWindowMinhashes.cend(); ++it)
          {
            auto offsetValues = it->second;

            if(offsetValues.first != NA)    //This hash occured in the reference window
              uniqueHashCount += 1;
          }

          float actualDensity = uniqueHashCount * 1.0 / len;
          float expectedDensity = 2.0 / this->windowSize;

          referenceDNAComplexity = actualDensity/expectedDensity;

          statValues.push_back(referenceDNAComplexity);
        }
      }
  };

}

#endif
