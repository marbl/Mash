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
#include "base_types.hpp"
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
    private:

      //algorithm parameters
      const skch::Parameters &param;

      //reference sketch
      const skch::Sketch &refSketch;

      //(seq. id, start ref. position, end position) for L1
      typedef std::tuple<seqno_t, offset_t, offset_t> candidateLocus_t;

      //Container type for saving read sketches during L1 and L2 both
      typedef Sketch::MI_Type MinVec_Type;

      typedef Sketch::MIIter_t MIIter_t;

      //  (seq. id, begin iterator over minimizerIndex, end iterator over minimizerIndex, 
      //  equal minhash count, mapping strand) for L2
      //  [begin iterator, end iterator) represents the mapped region on reference
      typedef std::tuple<seqno_t, MIIter_t, MIIter_t, int, strand_t> mapLocus_t;

      const offset_t NA = std::numeric_limits<offset_t>::max();    //hash 'Not Available' marker

      //Slide read within superwindow by these many minimizers
      int stepMinCount;

      typedef std::function< void(const MappingResult&) > PostProcessResultsFn_t;

      //Custom function for post processing the results 
      PostProcessResultsFn_t processMappingResults;

    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] refSketch   reference sketch
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      Map(const skch::Parameters &p, const skch::Sketch &refsketch,
          PostProcessResultsFn_t f = nullptr) :
        param(p),
        refSketch(refsketch),
        processMappingResults(f)
    {
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

        std::ofstream outstrm(param.outFileName);

        for(const auto &fileName : param.querySequences)
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
            //Get optimal window size level depending upon the read length
            wsize_t optimalWindowSizeLevel;

            //Is the read too short?
            if(len < param.baseWindowSize || len < param.kmerSize || len < param.minReadLength)
            {
              seqCounter++;

#ifdef DEBUG
              std::cout << "WARNING, skch::Map::mapQuery, read is not long enough for mapping" << std::endl;
#endif

              continue;
            }
            else 
            {
              QueryMetaData <decltype(seq), MinVec_Type> Q;

              Q.seq = seq;
              Q.seqCounter = seqCounter;
              Q.len = len; 

              //Map this sequence
              bool mappingReported = mapSingleQuerySeq(Q, outstrm);
              if(mappingReported)
                totalReadsMapped++;

              seqCounter++;
              totalReadsPickedForMapping++;
            }
          }

          //Close the input file
          kseq_destroy(seq);  
          gzclose(fp);  
        }

        std::cout << "INFO, skch::Map::mapQuery, [count of mapped reads, reads qualified for mapping, total input reads] = [" << totalReadsMapped << ", " << totalReadsPickedForMapping << ", " << seqCounter << "]" << std::endl;

      }

      /**
       * @brief               map the parsed query sequence (L1 and L2 mapping)
       * @param[in] Q         metadata about query sequence
       * @param[in] outstrm   outstream stream where mappings will be reported
       */
      template<typename Q_Info>
        inline bool mapSingleQuerySeq(Q_Info &Q, std::ofstream &outstrm)
        {
#if ENABLE_TIME_PROFILE_L1_L2
          auto t0 = skch::Time::now();
#endif

          //Compute optimal window size for sketching this read
          if(param.staticWin)
            Q.optimalWindowSizeLevel = 0;
          else
            Q.optimalWindowSizeLevel = Stat::recommendedWindowLevelForRead(param.p_value,
                param.kmerSize, param.alphabetSize,
                param.percentageIdentity,
                Q.len, param.referenceSize,
                param.baseWindowSize, param.dynamicWinLevels);

#if ENABLE_TIME_PROFILE_L1_L2
          auto t1 = skch::Time::now();
#endif
          //L1 Mapping
          std::vector<candidateLocus_t> l1Mappings; 
          doL1Mapping(Q, l1Mappings);

#if ENABLE_TIME_PROFILE_L1_L2
          std::chrono::duration<double> timeSpentL1 = skch::Time::now() - t1;
          t1 = skch::Time::now();
#endif

          //L2 Mapping
          MappingResultsVector_t l2Mappings;
          bool mappingReported = doL2Mapping(Q, l1Mappings, l2Mappings);

          //Write mapping results to file
          reportL2Mappings(l2Mappings, outstrm);


#if ENABLE_TIME_PROFILE_L1_L2
          {
            std::chrono::duration<double> timeSpentL2 = skch::Time::now() - t1;
            std::chrono::duration<double> timeSpentMappingRead = skch::Time::now() - t0;
            int countL1Candidates = l1Mappings.size();

            std::cerr << Q.seq->name.s << " " << Q.len
              << " " << countL1Candidates 
              << " " << timeSpentL1.count() 
              << " " << timeSpentL2.count()
              << " " << timeSpentMappingRead.count()
              << "\n";
          }
#endif

          return mappingReported;
        }

      /**
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is 
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an 
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details 
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename Vec>
        void doL1Mapping(Q_Info &Q, Vec &l1Mappings)
        {
          //Vector of positions of all the hits 
          std::vector<MinimizerMetaData> seedHitsL1;

          ///1. Compute the minimizers

          CommonFunc::addMinimizers(Q.minimizerTableQuery, Q.seq, param.kmerSize, param.baseWindowSize, Q.optimalWindowSizeLevel, param.alphabetSize);

#ifdef DEBUG
          std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", minimizer count = " << Q.minimizerTableQuery.size() << "\n";
#endif

          ///2. Find the hits in the reference, pick 's' unique minimizers as seeds, 

          std::sort(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::lessByHash);

          //note : unique preserves the original relative order of elements 
          auto uniqEndIter = std::unique(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::equalityByHash);

          //This is the sketch size for estimating jaccard
          Q.sketchSize = std::distance(Q.minimizerTableQuery.begin(), uniqEndIter);

          int totalMinimizersPicked = 0;

          for(auto it = Q.minimizerTableQuery.begin(); it != uniqEndIter; it++)
          {
            //Check if hash value exists in the reference lookup index
            auto seedFind = refSketch.minimizerPosLookupIndex.find(it->hash);

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

          //Remove hits in the reference with window size level < optimalWindowSizeLevel
          seedHitsL1.erase( std::remove_if( seedHitsL1.begin(), seedHitsL1.end(),
                [&](MinimizerMetaData &e){
                return (e.w_lev < Q.optimalWindowSizeLevel); }),
              seedHitsL1.end()); 

          int minimumHits = Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity);

          this->computeL1CandidateRegions(Q, seedHitsL1, minimumHits, l1Mappings);

#ifdef DEBUG
          std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", Count of L1 hits in the reference = " << seedHitsL1.size() << ", minimum hits required for a candidate = " << minimumHits << ", Count of L1 candidate regions = " << l1Mappings.size() << "\n";

          for(auto &e : l1Mappings)
            std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", L1 candidate : [" << this->refSketch.metadata[std::get<0>(e)].name << " " << this->refSketch.metadata[std::get<0>(e)].len << " " << std::get<1>(e) << " " << std::get<2>(e) << "]\n";
#endif

        }

      /**
       * @brief                     Helper function to doL1Mapping()
       * @param[in]   Q             query
       * @param[in]   seedHitsL1    minimizer hits in the reference
       * @param[in]   minimumHits   estimated minimum hits required for significant match
       * @param[out]  l1Mappings    all the read mapping locations
       */
      template <typename Q_Info, typename Vec1, typename Vec2>
        void computeL1CandidateRegions(Q_Info &Q, Vec1 &seedHitsL1, int minimumHits, Vec2 &l1Mappings)
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
              if(it2->seqId == it->seqId && it2->pos - it->pos < Q.len)
              {
                //Save <1st pos --- 2nd pos>
                candidateLocus_t candidate(it->seqId, 
                    std::max(0, it2->pos - Q.len), it->pos + Q.len);

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
       * @brief                                 Revise L1 candidate regions to more precise locations
       * @param[in]   Q                         query sequence information
       * @param[in]   l1Mappings                candidate regions for query sequence found at L1
       * @param[out]  l2Mappings                Mapping results in the L2 stage
       * @return      T/F                       true if atleast 1 mapping region is proposed
       */
      template <typename Q_Info, typename VecIn, typename VecOut>
        bool doL2Mapping(Q_Info &Q, VecIn &l1Mappings, VecOut &l2Mappings)
        {
          //Range of sketch in query
          auto uniqEndIter = std::next(Q.minimizerTableQuery.begin(), Q.sketchSize);

          /*
           * Type of container used for walking read over the reference superwindow
           * We choose a std map : 
           *    [hash value -> (location in the reference superwindow, offset in query)]
           *    offset is multiplied by strand type (+1,-1) to embed this information
           *
           * map preserves the sorted order of hashes, so computing jaccard similarity 
           * is a linear walk, each time for the sliding window
           */

          //Step size
          //TODO: Implement efficient sliding while considering all read sized windows
          stepMinCount = getStepMinCount();

          typedef std::map< hash_t, std::pair<offset_t, offset_t> > slidingMapType;

          slidingMapType slidingWindowMinhashes;

          for(auto it = Q.minimizerTableQuery.begin(); it != uniqEndIter; it++)
            slidingWindowMinhashes.emplace(it->hash, std::make_pair(NA, (it->pos)*(it->strand)));     //[hash value] -> (NA, offset in query)

          bool mappingReported = false;

          ///2. Walk the read over the candidate regions and compute the jaccard similarity with minimum s sketches
          for(auto &candidateLocus: l1Mappings)
          {
            //Reset the map to just have the query minimizers
            slidingMapType slidingWindowMinhashesCpy = slidingWindowMinhashes;

            mapLocus_t l2;
            computeL2MappedRegions(Q, slidingWindowMinhashesCpy, candidateLocus, l2);

            //Compute mash distance using calculated jaccard
            float mash_dist = Stat::j2md(1.0 * std::get<3>(l2)/Q.sketchSize, param.kmerSize);

            //Compute lower bound to mash distance within 90% confidence interval
            float mash_dist_lower_bound = Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, 0.9);

            float nucIdentity = 100 * (1 - mash_dist);
            float nucIdentityUpperBound = 100 * (1 - mash_dist_lower_bound);

            /*
             * Compute addtional statistics to filter out false mappings
             * 1. Unique Minimizers in the reference
             */
            std::vector<float> mappingStatistics;

            slidingWindowMinhashesCpy = slidingWindowMinhashes;
            computeMappingStatistics(Q, slidingWindowMinhashesCpy, 
                std::get<1>(l2),std::get<2>(l2), 
                mappingStatistics); 

            /*    An alignment is reported if 
             *    the nucleotide identity is >= the percentage identity threshold
             *    estimated complexity (uniqueness of minimizers) is above 75%
             */
            if(nucIdentityUpperBound >= param.percentageIdentity && mappingStatistics[0] >= 0.75)
            {
              MappingResult res;

              //Save the output
              {
                res.queryLen = Q.len;
                res.refStartPos = std::get<1>(l2)->pos ;
                res.refEndPos = std::get<1>(l2)->pos + Q.len - 1;
                res.refSeqId = std::get<0>(l2);
                res.querySeqId = Q.seqCounter;
                res.nucIdentity = nucIdentity;
                res.nucIdentityUpperBound = nucIdentityUpperBound;
                res.sketchSize = Q.sketchSize;
                res.conservedSketches = std::get<3>(l2);
                res.strand = std::get<4>(l2);
                res.mappedRegionComplexity = mappingStatistics[0];
                res.queryName = Q.seq->name.s; 

                l2Mappings.push_back(res);
              }

              mappingReported = true;
            }
          }

          return mappingReported;
        }

      /**
       * @brief       Helper function to doL2Mapping()
       * @param[in]   Q                         query sequence information
       * @param[in]   slidingWindowMinhashes    container for computing min s hashes along the sliding read 
       * @param[in]   candidateLocus            candidate region computed at L1 stage
       * @param[out]  l2Mapping                 best read mapping location inside the given candidateLocus
       */
      template <typename Q_Info, typename MapT>
        void computeL2MappedRegions(Q_Info &Q,
            MapT &slidingWindowMinhashes, 
            candidateLocus_t &candidateLocus, mapLocus_t &l2Mapping)
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
          auto lastPosition = std::make_pair(refSequenceId, currentStart + Q.len - 1);

#ifdef DEBUG
          if(currentStart + Q.len - 1 > std::get<2>(candidateLocus))
            std::cout << "WARNING, skch::Map::computeL2MappedRegions, candidate region shorter than read length" << std::endl;
#endif

          /*
           * [l_iter, u_iter) represents the superwindow iterator range of minimizers 
           * in the reference sketch index against which jaccard similarity of the 
           * query is estimated
           *
           * After each window shift, we advance the range [l_iter, u_iter) forward
           * depending on the 'stepMinCount'
           *
           * Therefore, each shift requires 
           *      1. deletion of [l_iter (previous), l_iter) elements from the window
           *      2. addition of [u_iter (previous), u_iter) elements into the window
           */

          auto iterRange = this->refSketch.getIndexRange(startPosition, lastPosition);
          auto l_iter = iterRange.first;
          auto u_iter = iterRange.second;

          assert(std::distance(l_iter, u_iter) > 0);

          //Save the best mapping locations here
          MIIter_t bestStartIter = l_iter, bestEndIter = u_iter;
          int maxSharedMinimizers = 0;
          strand_t bestMatchStrandType;

          //Previous range is null at the beginning
          auto l_prevIter = l_iter;
          auto u_prevIter = l_iter;

          //Keep sliding till query crosses the candidate's end position
          //Check if u_iter offset is less than L1 candidate boundary, and reference sequence id should match
          while(u_iter->pos < std::get<2>(candidateLocus) && 
              std::prev(u_iter)->seqId == std::get<0>(candidateLocus))
          {
            assert(std::distance(l_iter , u_iter) > 0);

            /// 1. Push minimizers in the [u_prevIter, u_iter) range to the window

            for(auto it = u_prevIter; it != u_iter; it++)
            {
              hash_t hashvalue = it->hash;
              offset_t offsetinReference = it->pos;
              wsize_t windowSizeLevel = it->w_lev;
              strand_t strandType = it->strand;

              //Minimizer in the reference only participates if it is sketched for window >= optimalWindowSizeLevel
              if(windowSizeLevel >= Q.optimalWindowSizeLevel)
              {
                //if hash doesn't exist in window
                if(slidingWindowMinhashes.find(hashvalue) == slidingWindowMinhashes.end())
                  slidingWindowMinhashes[hashvalue] = std::make_pair(offsetinReference * strandType, NA);   //add the hash to window
                else
                  slidingWindowMinhashes[hashvalue].first = offsetinReference * strandType;                 //just revise the offset 
              }
            }

            /// 2. Remove minimizers in the [l_prevIter, l_iter) range from the window 

            for(auto it = l_prevIter; it != l_iter; it++)
            {
              hash_t hashvalue = it->hash;
              offset_t offsetinReference = it->pos;
              wsize_t windowLevel = it->w_lev;
              strand_t strandType = it->strand;

              //Minimizer in the reference only participates if it is sketched for window >= optimalWindowSizeLevel
              if(windowLevel >= Q.optimalWindowSizeLevel)
              {
                //hash must exist in the map, by our logic
                assert(slidingWindowMinhashes.find(hashvalue) != slidingWindowMinhashes.end());

                //only if hash and offset match
                if(slidingWindowMinhashes[hashvalue].first == offsetinReference * strandType) 
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
            }

            assert(slidingWindowMinhashes.size() >= Q.sketchSize);

            /// 3. Compute jaccard similarity now

            int uniqueHashes = 0, currentSharedMinimizers = 0, strandVotes = 0;

            for (auto it = slidingWindowMinhashes.cbegin(); it != slidingWindowMinhashes.cend(); ++it)
            {
              auto offsetValues = it->second;

              if(offsetValues.first != NA && offsetValues.second != NA)
              {
                currentSharedMinimizers += 1;

                if( (offsetValues.first ^ offsetValues.second) >= 0)
                  strandVotes += 1;   //Both strand types match
                else
                  strandVotes -= 1;   //Both strand types mis-match
              }

              uniqueHashes++;
              if(uniqueHashes == Q.sketchSize)
                break;
            }

            //Is it the best superwindow so far?
            if(currentSharedMinimizers > maxSharedMinimizers)
            {
              bestStartIter = l_iter;
              bestEndIter = u_iter;
              maxSharedMinimizers = currentSharedMinimizers;

              if(strandVotes > 0)
                bestMatchStrandType = strnd::FWD;
              else
                bestMatchStrandType = strnd::REV;
            }

            /// 4. Advance the iterator range

            //Shift the window to right by stepMinCount (set to 20)

            //Backup the previous iterator positions 
            l_prevIter = l_iter; u_prevIter = u_iter;

            //Compute the new iterator positions
            if (!this->refSketch.advanceIter(l_iter, stepMinCount) || !this->refSketch.advanceIter(u_iter, stepMinCount))
              break;
          }

          //Save the result
          l2Mapping = std::make_tuple(refSequenceId, bestStartIter, bestEndIter, maxSharedMinimizers, bestMatchStrandType);
        }

      /**
       * @brief                     Report the final L2 mappings to output stream
       * @param[in]   l2Mappings    mapping results
       * @param[in]   outstrm       file output stream object
       */
      void reportL2Mappings(MappingResultsVector_t &l2Mappings, 
          std::ofstream &outstrm)
      {
        float bestNucIdentity = 0;

        //Scan through the mappings to check best identity mapping
        for(auto &e : l2Mappings)
        {
          if(e.nucIdentity > bestNucIdentity)
            bestNucIdentity = e.nucIdentity;
        }

        //Print the results
        for(auto &e : l2Mappings)
        {
          //Report top 1% mappings (unless reportAll flag is true, in which case we report all)
          if(param.reportAll == true || e.nucIdentity >= bestNucIdentity - 1.0)
          {
            outstrm << e.queryName 
              << " " << e.queryLen 
              << " " << "0"
              << " " << e.queryLen - 1 
              << " " << (e.strand == strnd::FWD ? "+" : "-") 
              << " " << this->refSketch.metadata[e.refSeqId].name
              << " " << this->refSketch.metadata[e.refSeqId].len
              << " " << e.refStartPos 
              << " " << e.refEndPos
              << " " << e.nucIdentity;

            //Print some additional statistics
            outstrm << " " << e.conservedSketches 
              << " " << e.sketchSize 
              << " " << e.nucIdentityUpperBound
              << " " << e.mappedRegionComplexity;

            outstrm << "\n";

            //User defined processing of the results
            if(processMappingResults != nullptr)
              processMappingResults(e);
          }
        }
      }

      /**
       * @brief                               compute additional statistics for mapping (post-L2)
       * @param[in]   Q                       query sequence information
       * @param[in]   slidingWindowMinhashes    
       * @param[in]   refSequenceId           begin iterator on reference index by L2 mapping
       * @param[in]   refEndItr               end iterator on reference index by L2 mapping 
       * @param[out]  statValues     
       */
      template <typename Q_Info, typename MapT, typename statVec_t>
        void computeMappingStatistics(Q_Info &Q, 
            MapT &slidingWindowMinhashes, 
            MIIter_t refStartItr, MIIter_t refEndItr,
            statVec_t &statValues)
        {
          //Push minimizers in the [refStartItr, refEndItr) range to the map
          for(auto it = refStartItr; it != refEndItr; it++)
          {
            hash_t hashvalue = it->hash;
            offset_t offsetinReference = it->pos;

            //if hash doesn't exist in window
            if(slidingWindowMinhashes.find(hashvalue) == slidingWindowMinhashes.end())
              slidingWindowMinhashes[hashvalue] = std::make_pair(offsetinReference, NA);   //add the hash to window
            else
              slidingWindowMinhashes[hashvalue].first = offsetinReference;                 //just revise the offset 
          }

          ///1. Unique Minimizers in the reference (complexity)
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

            float actualDensity = uniqueHashCount * 1.0 / Q.len;
            float expectedDensity = 2.0 / param.baseWindowSize;

            referenceDNAComplexity = actualDensity/expectedDensity;

            statValues.push_back(referenceDNAComplexity);
          }
        }

      /**
       * @brief   sliding step size during L2 stage
       * @detail  Ideally we should evaluate all sliding windows by step 1, 
       *          which is theoretically feasible using STL map. But due to 
       *          multiple corner cases, we resort to do sliding the read by 
       *          more than 1 step.
       * @return  count of minimizers to jump while sliding 
       */
      inline int getStepMinCount()
      {
        float sketchDensity = 2.0/param.baseWindowSize;
        float percentReadSkip = 0.2;      //Jump 20% minimum read size
        return (param.minReadLength * percentReadSkip) * sketchDensity;
      }

    public:

      /**
       * @brief     An optional utility function to save the 
       *            reported results
       */
      static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
      {
        v.push_back(reportedL2Result);
      }
  };

}

#endif
