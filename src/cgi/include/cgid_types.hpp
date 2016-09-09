/**
 * @file    cgid_types.hpp
 * @brief   Critical type defintions for ANI algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CGID_TYPES_HPP 
#define CGID_TYPES_HPP

#include <tuple>
#include "base_types.hpp"

namespace cgi
{
  //Final mapping result
  struct MappingResult_CGI
  {
    skch::seqno_t refSeqId;           //internal sequence id of the reference contig
    skch::seqno_t querySeqId;         //name of query sequence
    float nucIdentity;                //calculated identity

    //Lexographical less than comparison
    bool operator <(const MappingResult_CGI& x) const {
      return std::tie(refSeqId, querySeqId, nucIdentity) 
        < std::tie(x.refSeqId, x.querySeqId, x.nucIdentity);
    }
  };

  struct CGI_Results
  {
    skch::seqno_t refSeqId;
    skch::seqno_t countSeq;
    float identity;

    //Default comparison is by identity
    bool operator <(const CGI_Results& x) const {
      return identity < x.identity;
    }
  };

}

#endif
