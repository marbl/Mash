/**
 * @file    base_types.hpp
 * @brief   Critical type defintions for mapping algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef BASE_TYPES_MAP_HPP 
#define BASE_TYPES_MAP_HPP

#include <tuple>

namespace skch
{
  typedef uint32_t hash_t;  //hash type
  typedef int offset_t;     //position within sequence
  typedef int seqno_t;      //sequence counter in file
  typedef int wsize_t;      //window size type 

  //C++ timer
  typedef std::chrono::high_resolution_clock Time;

  //Information about each minimizer
  struct MinimizerInfo
  {
    hash_t hash;            //hash value
    seqno_t seqId;          //sequence or contig id
    offset_t pos;           //position within sequence
    wsize_t win;            //associated window size

    bool operator <(const MinimizerInfo& x) {
      return std::tie(hash, seqId, pos, win) 
        < std::tie(x.hash, x.seqId, x.pos, x.win);
    }

    bool operator ==(const MinimizerInfo& x) {
      return std::tie(hash, seqId, pos, win) 
        == std::tie(x.hash, x.seqId, x.pos, x.win);
    }

    bool operator !=(const MinimizerInfo& x) {
      return std::tie(hash, seqId, pos, win) 
        != std::tie(x.hash, x.seqId, x.pos, x.win);
    }

    static bool equalityByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash == y.hash;
    }

    static bool lessByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash < y.hash;
    }

  };

  //Type for map value type used for
  //L1 stage lookup index
  struct MinimizerMetaData
  {
    seqno_t seqId;          //sequence or contig id
    offset_t pos;           //position within sequence
    wsize_t win;            //associated window size

    bool operator <(const MinimizerMetaData& x) const {
      return std::tie(seqId, pos, win) 
        < std::tie(x.seqId, x.pos, x.win);
    }
  };

  typedef hash_t MinimizerMapKeyType;
  typedef std::vector<MinimizerMetaData> MinimizerMapValueType;

  //Metadata recording for contigs in the reference DB
  struct ContigInfo
  {
    std::string name;       //Name of the sequence
    offset_t len;           //Length of the sequence
  };
}

#endif
