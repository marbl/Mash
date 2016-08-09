/**
 * @file    commonFunc.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMMON_FUNC_HPP 
#define COMMON_FUNC_HPP

#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>

#ifdef USE_BOOST
    #include <boost/math/distributions/binomial.hpp>
    using namespace::boost::math;
#else
    #include <gsl/gsl_cdf.h>
#endif

//Own includes
#include "map_parameters.hpp"
#include "winSketch.hpp"

//External includes
#include "murmur3.h"
#include "kseq.h"
#include "prettyprint.hpp"

namespace skch
{
  /**
   * @class     skch::CommonFunc
   * @brief     Implements frequently used common functions
   */
  namespace CommonFunc
  {
    //seed for murmerhash
    const int seed = 42;

    /**
     * @brief   reverse complement of kmer (borrowed from mash)
     */
    inline void reverseComplement(const char * src, char * dest, int length) 
    {
      for ( int i = 0; i < length; i++ )
      {    
        char base = src[i];

        switch ( base )
        {    
          case 'A': base = 'T'; break;
          case 'C': base = 'G'; break;
          case 'G': base = 'C'; break;
          case 'T': base = 'A'; break;
          default: break;
        }    

        dest[length - i - 1] = base;
      }    
    }

    template <typename KSEQ>
      inline void makeUpperCase(KSEQ *seq)
      {
        for ( int i = 0; i < seq->seq.l; i++ )
        {
          if (seq->seq.s[i] > 96 && seq->seq.s[i] < 123)
          {
            seq->seq.s[i] -= 32;
          }
        }
      }

    /**
     * @brief   hashing kmer string (borrowed from mash)
     */
    inline hash_t getHash(const char * seq, int length)
    {
      char data[16];
      MurmurHash3_x64_128(seq, length, seed, data);

      hash_t hash;

      hash = *((hash_t *)data);

      return hash;
    }

    /**
     * @brief       compute winnowed minimizers from a given sequence and add to the index
     * @param[out]  minimizerIndex  minimizer table storing minimizers and their position as we compute them
     * @param[in]   seq             kseq fasta/q parser
     * @param[in]   kmerSize
     * @param[in]   windowSize
     * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
     */
    template <typename T, typename KSEQ>
      inline void addMinimizers(std::vector<T> &minimizerIndex, KSEQ *seq, int kmerSize, int windowSize, seqno_t seqCounter)
      {
        //Double-ended queue (saves minimum at front end)
        std::deque< std::pair<hash_t, offset_t> > Q;

        makeUpperCase(seq);

        //length of the sequencd
        offset_t len = seq->seq.l;

        //Compute reverse complement of seq
        char *seqRev = new char[len];
        CommonFunc::reverseComplement(seq->seq.s, seqRev, len);


        for(offset_t i = 0; i < len - kmerSize + 1; i++)
        {
          //Hash kmers
          hash_t hashFwd = CommonFunc::getHash(seq->seq.s + i, kmerSize); 
          hash_t hashBwd = CommonFunc::getHash(seqRev + len - i - kmerSize, kmerSize);

          //Consider non-symmetric kmers only
          if(hashBwd != hashFwd)
          {
            //Take minimum value of kmer and its reverse complement
            hash_t currentKmer = std::min(hashFwd, hashBwd);

            //If front minimum is not in the current window, remove it
            while(!Q.empty() && Q.front().second <=  i - windowSize)
              Q.pop_front();

            //Hashes less than equal to currentKmer are not required
            //Remove them from Q (back)
            while(!Q.empty() && currentKmer <= Q.back().first) 
              Q.pop_back();

            //Push currentKmer into back of the queue
            Q.emplace_back(currentKmer, i); 

            //Select the minimizer from Q and put into index
            //Ignore the minimizers from first few incomplete sliding windows
            if(i >= windowSize - 1)
            {
              auto potentialMinimizer = std::make_tuple(Q.front().first, seqCounter, Q.front().second);

              //We save the minimizer if we are seeing it for first time
              if(minimizerIndex.empty() || minimizerIndex.back() != potentialMinimizer)
              {
                minimizerIndex.push_back(potentialMinimizer);
              }
            }
          }
        }

#ifdef DEBUG
        //std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted following minimizers: ";
        //for(auto &e: minimizerIndex) std::cout << "\n" << e;
        //std::cout << std::endl;
#endif

        delete [] seqRev;
      }

    /**
     * @brief       overloaded function for case where seq. counter does not matter
     */
    template <typename T, typename KSEQ>
      inline void addMinimizers(std::vector<T> &minimizerIndex, KSEQ *seq, int kmerSize, int windowSize)
      {
        addMinimizers(minimizerIndex, seq, kmerSize, windowSize, 0);
      }

    /**
     * @brief       Estimate minimum number of shared minimizers/kmers to achieve the desired identity
     */
    inline int estimateMinimumHits(int s, int kmerSize, int windowSize, float perc_identity)
    {
      //function to convert mash distance to jaccard similarity
      auto mashToJaccard = [&](float mash_dist){
        return 1.0 / ( 2.0 * exp(kmerSize * mash_dist) - 1.0);
      }; 

      //function to convert jaccard to min hits
      //Atleast these many minimizers should match for jaccard identity
      auto getMinimumSharedMinimizers = [&](float jaccard){
        return floor (1.0 * s * jaccard); 
      };

      //Compute the estimate
      float mash_dist = 1.0 - perc_identity/100.0;
      float jaccardSimilarity = mashToJaccard(mash_dist);

      return getMinimumSharedMinimizers(jaccardSimilarity);

    }

    /**
     * @brief                     Compute mash nucleotide identity given k, s, w
     * @param[in]   denom         sketch size
     * @param[in]   common        matched sketches
     * @param[in]   kmerSize      kmer size
     * @param[out]  nucIdentity   mash nucleotide identity [0-100]
     */
    inline void computeMashNucIdentity(int denom, int common, int kmerSize, float &nucIdentity)
    {
      assert(denom > 0);

      float jaccard = 1.0 * common / denom;

      if ( common == denom )
        nucIdentity = 100.0;

      else if( common == 0 )
        nucIdentity= 0.0;

      else
        nucIdentity = 100 * (1 + log(2.0 * jaccard / (1.0 + jaccard)) / kmerSize);
    }

    /**
     * @brief                             Compute upper bound on the identity (eq. to lower bound on mash distance)
     * @details                           Used to get a conservative identity estimate for filtering out poor L2 hits
     * @param[in]   identity              mash identity [0-100]
     * @param[in]   s                     sketch size
     * @param[in]   kmerSize              kmer size
     * @param[in]   prob                  probability that identity has an upper bound of upperBoundIdentity
     *                                    Higher probability implies larger bound
     * @param[out]  upperBoundIdentity    Upper bound [0 - 100]
     */
    inline void computeUpperBoundIdentity(float identity, int s, int kmerSize, float prob, float &upperBoundIdentity)
    {
      //function to convert mash distance to jaccard similarity
      auto mashToJaccard = [&](float mash_dist){
        return 1.0 / ( 2.0 * exp(kmerSize * mash_dist) - 1.0);
      }; 

      float q2 = (1.0 - prob)/2;

      float mash_dist = 1.0 - identity/100;

#ifdef USE_BOOST
      int x = quantile(complement(binomial(s, mashToJaccard(mash_dist)), q2));

#else   //GSL

      int x = s;
			while ( x > 0 )
      {
        double cdf_complement = gsl_cdf_binomial_Q(x-1, mashToJaccard(mash_dist), s);

				if ( cdf_complement > q2 )
				{
					break;
				}
				
				x--;
      }
#endif

      float je = float(x) / s;

      float distance =  (-1.0 / kmerSize) * log(2.0 * je / (1.0 + je)) ;

      upperBoundIdentity = 100 * (1-distance);
    }

    /**
     * @brief           Functor for comparing tuples by single index layer
     * @tparam layer    Tuple's index which is used for comparison
     * @tparam op       comparator, default as std::less
     */
    template <size_t layer, template<typename> class op = std::less>
      struct TpleComp
      {
        //Compare two tuples using their values
        template<typename T>
          bool operator() (T const &t1, T const &t2)
          {
            return op<typename std::tuple_element<layer, T>::type>() (std::get<layer>(t1), std::get<layer>(t2));
          }
      };

  }
}

#endif
