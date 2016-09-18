/**
 * @file    computeCoreIdentity.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CGI_IDENTITY_HPP 
#define CGI_IDENTITY_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>

//Own includes
#include "base_types.hpp"
#include "cgid_types.hpp"

//External includes
#include "prettyprint.hpp"

namespace cgi
{
  void computeCGI(skch::MappingResultsVector_t &results,
     skch::Sketch &refSketch,
     std::string &fileName)
  {
    std::vector<MappingResult_CGI> shortResults;
    shortResults.reserve(results.size());

    for(auto &e : results)
    {
      shortResults.push_back(MappingResult_CGI{
          e.refSeqId,
          e.querySeqId,
          e.nucIdentity
          });
    }

    std::sort(shortResults.begin(), shortResults.end());

    std::vector<MappingResult_CGI> singleQueryRefResults;

    for(auto &e : shortResults)
    {
      if(singleQueryRefResults.empty())
        singleQueryRefResults.push_back(e);

      else if ( !(
            e.refSeqId == singleQueryRefResults.back().refSeqId && 
            e.querySeqId == singleQueryRefResults.back().querySeqId))
        singleQueryRefResults.push_back(e);

      else
        singleQueryRefResults.back().nucIdentity = e.nucIdentity;
    }

    std::vector<cgi::CGI_Results> CGI_ResultsVector;


    for(auto it = singleQueryRefResults.begin(); it != singleQueryRefResults.end();)
    {
      skch::seqno_t currentRefId = it->refSeqId;

      auto rangeEndIter = std::find_if(it, singleQueryRefResults.end(), [&](const MappingResult_CGI& e) 
          { 
            return e.refSeqId != currentRefId; 
          } );

      float sumIdentity = 0.0;

      for(auto it2 = it; it2 != rangeEndIter; it2++)
      {
        sumIdentity += it2->nucIdentity;
      }

      CGI_Results currentResult;
      currentResult.refSeqId = currentRefId;
      currentResult.countSeq = std::distance(it, rangeEndIter);
      currentResult.identity = sumIdentity/currentResult.countSeq;

      CGI_ResultsVector.push_back(currentResult);

      //Advance the iterator it
      it = rangeEndIter;
    }

    //Sort in decreasing order of matches
    std::sort(CGI_ResultsVector.rbegin(), CGI_ResultsVector.rend());

    std::ofstream outstrm(fileName);

    for(auto &e : CGI_ResultsVector)
    {
      outstrm << refSketch.metadata[e.refSeqId].name
        << " " << e.identity 
        << " " << e.countSeq
        << " " << refSketch.metadata[e.refSeqId].len
        << "\n";
    }
  }

}

#endif
