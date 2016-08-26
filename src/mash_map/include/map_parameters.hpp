/**
 * @file    map_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_CONFIG_HPP 
#define SKETCH_CONFIG_HPP

#include <vector>

//Own includes

namespace skch
{
  typedef uint32_t hash_t;  //hash type
  typedef int offset_t;     //position within sequence
  typedef int seqno_t;      //sequence counter in file
  typedef int wsize_t;      //window size type 

  //C++ timer
  typedef std::chrono::high_resolution_clock Time;

  /**
   * @brief   configuration parameters for building sketch
   *          expected to be initialized using command line arguments
   */
  struct Parameters
  {
    int kmerSize;                                     //kmer size for sketching

    int baseWindowSize;                               //base window size used for sketching 
    bool dynamicWin;                                  //switch for dynamic windowing
    int dynamicWinLevels;                             //no. of hierarchical levels for dynamic windowing
                                                      // including base level

    int minReadLength;                                //minimum read length which code maps
    int alphabetSize;                                 //alphabet size
    uint64_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    double p_value;                                   //user defined threshold for p value
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
  };
}

#endif
