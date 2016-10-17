/**
 * @file    winSketch.hpp
 * @brief   routines to index the reference 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WIN_SKETCH_HPP 
#define WIN_SKETCH_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <cassert>
#include <zlib.h>  

//Own includes
#include "base_types.hpp"
#include "map_parameters.hpp"
#include "commonFunc.hpp"

//External includes
#include "kseq.h"
#include "murmur3.h"
#include "prettyprint.hpp"

KSEQ_INIT(gzFile, gzread)

namespace skch
{
  /**
   * @class     skch::Sketch
   * @brief     sketches and indexes the reference (subject sequence)
   * @details  
   *            1.  Minimizers are computed in streaming fashion
   *                Computing minimizers is using double ended queue which gives
   *                O(reference size) complexity
   *                Algorithm described here:
   *                https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
   *
   *            2.  Index hashes into appropriate format to enable fast search at L1 mapping stage
   */
    class Sketch
    {
      //private members
    
      //algorithm parameters
      const skch::Parameters &param;

      //Ignore top % most frequent minimizers while lookups
      const float percentageThreshold = 0.001;

      //Minimizers that occur this or more times will be ignored (computed based on percentageThreshold)
      int freqThreshold = std::numeric_limits<int>::max();

      //Make the default constructor private, non-accessible
      Sketch();

      public:

      typedef std::vector< MinimizerInfo > MI_Type;
      using MIIter_t = MI_Type::const_iterator;

      //Keep sequence length, name that appear in the sequence (for printing the mappings later)
      std::vector< ContigInfo > metadata;

      /*
       * Keep the information of what sequences come from what file#
       * Example [a, b, c] implies 
       *  file 0 contains 0 .. a-1 sequences
       *  file 1 contains a .. b-1 
       *  file 2 contains b .. c-1
       */
      std::vector< int > sequencesByFileInfo;

      //Index for fast seed lookup
      /*
       * [minimizer #1] -> [pos1, pos2, pos3 ...]
       * [minimizer #2] -> [pos1, pos2...]
       * ...
       */
      using MI_Map_t = std::unordered_map< MinimizerMapKeyType, MinimizerMapValueType >;
      MI_Map_t minimizerPosLookupIndex;

      private:

      /**
       * Keep list of minimizers, sequence# , their position within seq , here while parsing sequence 
       * Note : position is local within each contig
       * Hashes saved here are non-unique, ordered as they appear in the reference
       */
      MI_Type minimizerIndex;

      //Frequency histogram of minimizers
      //[... ,x -> y, ...] implies y number of minimizers occur x times
      std::map<int, int> minimizerFreqHistogram;

      public:

      /**
       * @brief   constructor
       *          also builds, indexes the minimizer table
       */
      Sketch(const skch::Parameters &p) 
        :
          param(p) {
            this->build();
            this->index();
            this->computeFreqHist();
          }

      private:

      /**
       * @brief     build the sketch table
       * @details   compute and save minimizers from the reference sequence(s)
       *            assuming a fixed window size
       */
      void build()
      {

        //sequence counter while parsing file
        seqno_t seqCounter = 0;

        for(const auto &fileName : param.refSequences)
        {

#ifdef DEBUG
        std::cout << "INFO, skch::Sketch::build, building minimizer index for " << fileName << std::endl;
#endif

          //Open the file using kseq
          FILE *file = fopen(fileName.c_str(), "r");
          gzFile fp = gzdopen(fileno(file), "r");
          kseq_t *seq = kseq_init(fp);


          //size of sequence
          offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {
            //Save the sequence name
            metadata.push_back( ContigInfo{seq->name.s, (offset_t)seq->seq.l} );

            //Is the sequence too short?
            if(len < param.baseWindowSize || len < param.kmerSize)
            {
#ifdef DEBUG
              std::cout << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer and window size" << std::endl;
#endif
              seqCounter++;
              continue;  
            }
            else
            {
              wsize_t baseWindowLevel = 0;  //Perform winnowing with baseWindowSize

              skch::CommonFunc::addMinimizers(this->minimizerIndex, seq, param.kmerSize, param.baseWindowSize, baseWindowLevel, param.alphabetSize, seqCounter);
            }

            seqCounter++;
          }

          sequencesByFileInfo.push_back(seqCounter);

          kseq_destroy(seq);  
          gzclose(fp); //close the file handler 
          fclose(file);
        }

        std::cout << "INFO, skch::Sketch::build, minimizers picked from reference = " << minimizerIndex.size() << std::endl;

        if(param.staticWin == false)
          reviseTableDynamicWindows();

      }

      /**
       * @brief   revise the sketch table
       */
      void reviseTableDynamicWindows()
      {
        int baseWindowSize = param.baseWindowSize;
        int maxWindowSize = pow(2, param.dynamicWinLevels - 1) * baseWindowSize;

        //Revise for window sizes 2^0 * baseWindowSize ... 2^(dynamicWinLevels - 1) * baseWindowSize
        for(int winLevel = 1; winLevel <= param.dynamicWinLevels - 1; winLevel++)
        {
          //Dynamic winnowing for window size 'win'
          int w = pow(2, winLevel) * param.baseWindowSize;

          // Double-ended queue (saves minimum minimizer in 'win' sized window
          // at front end)
          std::deque< MI_Type::iterator > Q;

          for(auto it = this->minimizerIndex.begin(); it != minimizerIndex.end(); )
          {
            //Process needs to be done in chunks for each sequence
            //Because sliding windows should not overlap across 2 sequences
            auto currentSequenceRangeStart = it;
            seqno_t currentSequenceId = currentSequenceRangeStart->seqId; 

            auto currentSequenceRangeEnd = currentSequenceRangeStart;

            //Find end of current sequence block
            while( currentSequenceRangeEnd != minimizerIndex.end() && 
                currentSequenceRangeEnd->seqId == currentSequenceId)
            {
              currentSequenceRangeEnd++;
            }

            /*
             * [currentSequenceRangeStart, currentSequenceRangeEnd) represents current 
             * sequence range in the minimizer index
             */
            for(auto it2 = currentSequenceRangeStart; it2 != currentSequenceRangeEnd; it2++)
            {
              offset_t pos = it2->pos;
              hash_t hash_val = it2->hash;

              //If front minimum is not in the current window, remove it
              while(!Q.empty() && Q.front()->pos <=  pos - w)
                Q.pop_front();

              //  Minimizers less than equal to current minimizer 
              //  are not required.
              //  Remove them from Q (back)
              while(!Q.empty() && Q.back()->hash >= hash_val) 
                Q.pop_back();

              //Push current minimizer into back of the queue
              Q.emplace_back(it2); 

              if(pos >= w - 1)
              {
                //Pick the minimum minimizer
                //Update its window level
                Q.front()->w_lev = winLevel;
              }
            
            }

            Q.clear();

            it = currentSequenceRangeEnd;
            
          }
        }

        std::cout << "INFO, skch::Sketch::reviseTableDynamicWindows, updated minimizers for dynamic windowing : w = [" << baseWindowSize <<  " ... " << maxWindowSize << "]" << std::endl;

#ifdef DEBUG
        for(int winLevel = 1; winLevel <= param.dynamicWinLevels -1; winLevel++)
        {
          uint64_t countMinw = 0;
          for(auto &e : minimizerIndex)
            if(e.w_lev >= winLevel)
              countMinw++;
          std::cout << "Count of minimizers associated with window size " << winLevel * param.baseWindowSize << " = " << countMinw << std::endl; 
        }
#endif
      }

      /**
       * @brief   build the index for fast lookups using minimizer table
       */
      void index()
      {
        //Parse all the minimizers and push into the map
        for(auto &e : minimizerIndex)
        {
          // [hash value -> info about minimizer]
          minimizerPosLookupIndex[e.hash].push_back( 
              MinimizerMetaData{e.seqId, e.pos, e.w_lev, e.strand});
        }

        std::cout << "INFO, skch::Sketch::index, unique minimizers = " << minimizerPosLookupIndex.size() << std::endl;
      }

      /**
       * @brief   report the frequency histogram of minimizers using position lookup index
       *          and compute which high frequency minimizers to ignore
       */
      void computeFreqHist()
      {

        //1. Compute histogram

        for(auto &e : this->minimizerPosLookupIndex)
          this->minimizerFreqHistogram[e.second.size()] += 1;

        std::cout << "INFO, skch::Sketch::computeFreqHist, Frequency histogram of minimizers = " <<  *this->minimizerFreqHistogram.begin() <<  " ... " << *this->minimizerFreqHistogram.rbegin() << std::endl;

        //2. Compute frequency threshold to ignore most frequent minimizers

        int64_t totalUniqueMinimizers = this->minimizerPosLookupIndex.size();
        int64_t minimizerToIgnore = totalUniqueMinimizers * percentageThreshold / 100;

        int64_t sum = 0;

        //Iterate from highest frequent minimizers
        for(auto it = this->minimizerFreqHistogram.rbegin(); it != this->minimizerFreqHistogram.rend(); it++)
        {
          sum += it->second; //add frequency
          if(sum < minimizerToIgnore)
          {
            this->freqThreshold = it->first;
            //continue
          }
          else if(sum == minimizerToIgnore)
          {
            this->freqThreshold = it->first;
            break;
          }
          else
          {
            break;
          }
        }

        if(this->freqThreshold != std::numeric_limits<int>::max())
          std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "\%, ignore minimizers occurring >= " << this->freqThreshold << " times during lookup." << std::endl;
        else
          std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "\%, consider all minimizers during lookup." << std::endl;

      }

      public:

      int getFreqThreshold() const
      {
        return this->freqThreshold;
      }

      /*
       * @brief         get [begin, end) range on minimizerIndex given
       *                start and last positions of superwindow 
       * @param[in]     start     start position of superwindow in the reference
       * @param[in]     last      last position of superwindow in the reference
       * @return        pair of const iterators [begin, end) range on minimizerIndex
       *                enclosing the start and last positions
       */
      template<class POS>
        std::pair<MIIter_t, MIIter_t> getIndexRange(POS start, POS last) const
        {
          //Sanity check for type POS
          static_assert(std::is_same<POS, std::pair<seqno_t, offset_t> >::value, "Type check failed");
          assert(last > start);

          MIIter_t begin = std::lower_bound(this->minimizerIndex.begin(), this->minimizerIndex.end(), start, cmp);
          MIIter_t end = std::upper_bound(this->minimizerIndex.begin(), this->minimizerIndex.end(), last, cmp);

          return std::make_pair(begin, end);
        }

      /*
       * @brief                   count unique minimizer values in the given range
       *                          on minimizerIndex
       * @param[in]     start     start position of superwindow in the reference
       * @param[in]     last      last position of superwindow in the reference
       * @return                  unique count 
       */
      template<class POS>
        uint64_t countUniqueMinimizers(POS start, POS last) const
        {
          //Sanity check for type POS
          static_assert(std::is_same<POS, std::pair<seqno_t, offset_t> >::value, "Type check failed");
          assert(last > start);

          std::vector<hash_t> allHashes;

          //Compute the begin and end iterators corresponding to given start, last positions
          MIIter_t begin = std::lower_bound(this->minimizerIndex.begin(), this->minimizerIndex.end(), start, cmp);
          MIIter_t end = std::upper_bound(this->minimizerIndex.begin(), this->minimizerIndex.end(), last, cmp);

          //Parse through all the minimizers
          for(auto it = begin; it != end; it++)
            allHashes.push_back(it->hash);

          //Compute the count of unique minimizers
          std::sort(allHashes.begin(), allHashes.end());
          auto unique_end = std::unique(allHashes.begin(), allHashes.end());

          return std::distance(allHashes.begin(), unique_end);
        }

      /*
       * @brief                   advance iterator over minimizerIndex 
       *                          clamp to end if shift is off limits
       * @param[in] itr           iterator over minimizerIndex
       * @param[in] n             offset by which to shift itr
       * @return                  true if itr is not clamped to end. So true also implies that 
       *                          itr can be dereferenced safely
       */
      inline bool advanceIter(MIIter_t &itr, std::size_t n) const
      {
        std::size_t maxShift = std::distance(itr, this->minimizerIndex.end());

        if(n < maxShift)
        {
          itr =  std::next(itr, n);
          return true;
        }
        else
        {
          itr = this->minimizerIndex.end();
          return false;
        }
      }



      private:

      /*
       * @brief     functor for comparing minimizers by their position in minimizerIndex
       * @details   used for locating minimizers using known position range in the reference
       *            helper function to getIndexRange
       */
      struct compareMinimizersByPos
      {
        typedef std::pair<seqno_t, offset_t> P;

        bool operator() (const MinimizerInfo &m, const P &val)
        {
          P minimizerPosition(m.seqId, m.pos);
          return (minimizerPosition < val);
        }

        bool operator() (const P &val, const MinimizerInfo &m)
        {
          P minimizerPosition(m.seqId, m.pos);
          return (val < minimizerPosition);
        }
      } cmp;

    };
}

#endif
