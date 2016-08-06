/**
 * @file    parseCmdArgs.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PARSE_CMD_HPP 
#define PARSE_CMD_HPP

#include <iostream>
#include <string>
#include <fstream>

//Own includes
#include "map_parameters.hpp"

//External includes
#include "argvparser.hpp"

namespace skch
{
  /*
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd)
  {
    cmd.setIntroductoryDescription("Approximate read mapper based on Jaccard similarity");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("subject", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subject","s");

    cmd.defineOption("subjectList", "a file containing list of reference files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subjectList","sl");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("queryList", "a file containing list of query files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("queryList","ql");

    cmd.defineOption("kmer", "kmer size <= 16 [default 16]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("window", "window size [default : 100]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("window","w");

    cmd.defineOption("minReadLen", "minimum read length to map [default : 5000]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("minReadLen","m");

    cmd.defineOption("perc_identity", "threshold for identity [default : 85]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("perc_identity","pi");

    cmd.defineOption("output", "output file name", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");
  }

  /*
   * @brief                   Parse the file containing path to reference or query files
   * @param[in]   fileToRead
   * @param[out]  fileList    
   */
  template <typename VEC>
    void parseFileList(std::string &fileToRead, VEC &fileList)
    {
      std::string line;

      std::ifstream in(fileToRead);

      if (in.fail())
      {
        std::cerr << "ERROR, skch::parseFileList, Could not open " << fileToRead << "\n";
        exit(1);
      }

      while (std::getline(in, line))
      {
        fileList.push_back(line);
      }
    }

  /*
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      CommandLineProcessing::ArgvParser &cmd, 
      skch::Parameters &parameters)
  {
    int result = cmd.parse(argc, argv);

    //Make sure we get the right command line args
    if (result != ArgvParser::NoParserError)
    {
      std::cout << cmd.parseErrorDescription(result) << "\n";
      exit(1);
    }
    else if (!cmd.foundOption("subject") && !cmd.foundOption("subjectList"))
    {
      std::cout << "Provide reference file (s)\n";
      exit(1);
    }
    else if (!cmd.foundOption("query") && !cmd.foundOption("queryList"))
    {
      std::cout << "Provide reference file (s)\n";
      exit(1);
    }

    std::stringstream str;

    //Parse reference files
    if(cmd.foundOption("subject"))
    {
      std::string ref;

      str << cmd.optionValue("subject");
      str >> ref;

      parameters.refSequences.push_back(ref);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("subjectList");
      str >> listFile;

      parseFileList(listFile, parameters.refSequences);
    }

    str.clear();

    //Parse query files
    if(cmd.foundOption("query"))
    {
      std::string query;

      str << cmd.optionValue("query");
      str >> query;

      parameters.querySequences.push_back(query);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("queryList");
      str >> listFile;

      parseFileList(listFile, parameters.querySequences);
    }

    str.clear();

    //Parse algorithm parameters
    if(cmd.foundOption("kmer"))
    {
      str << cmd.optionValue("kmer");
      str >> parameters.kmerSize;
      str.clear();
    }
    else
      parameters.kmerSize = 16;

    if(cmd.foundOption("window"))
    {
      str << cmd.optionValue("window");
      str >> parameters.windowSize;
      str.clear();
    }
    else
      parameters.windowSize = 100;

    if(cmd.foundOption("minReadLen"))
    {
      str << cmd.optionValue("minReadLen");
      str >> parameters.minReadLength;
      str.clear();
    }
    else
      parameters.minReadLength = 5000;

    if(cmd.foundOption("perc_identity"))
    {
      str << cmd.optionValue("perc_identity");
      str >> parameters.percentageIdentity;
      str.clear();
    }
    else
      parameters.percentageIdentity = 85;

    str << cmd.optionValue("output");
    str >> parameters.outFileName;
    str.clear();

#ifdef DEBUG
    std::cout << "Reference = " << parameters.refSequence << std::endl;
    std::cout << "Query = " << parameters.querySequence << std::endl;
    std::cout << "Kmer size = " << parameters.kmerSize << std::endl;
    std::cout << "Window size = " << parameters.windowSize << std::endl;
    std::cout << "Percentage identity threshold = " << parameters.percentageIdentity << std::endl;
    std::cout << "Mapping output saved to file " << parameters.outFileName << std::endl;
#endif

  }
}


#endif
