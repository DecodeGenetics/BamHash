#include <iostream>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <seqan/arg_parse.h>

#include "bamhash_checksum_common.h"

struct Fastainfo {
  std::vector<std::string> fastafiles;
  bool debug;
  bool noReadNames;

  Fastainfo() : debug(false), noReadNames(false) {}

};

seqan::ArgumentParser::ParseResult
parseCommandLine(Fastainfo& options, int argc, char const **argv) {
  // Setup ArgumentParser.
  seqan::ArgumentParser parser("bamhash_checksum_fasta");

  setShortDescription(parser, "Checksum of a set of fasta files");
  setVersion(parser, BAMHASH_VERSION);
  setDate(parser, "Okt 2018");

  addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI<in1.fasta>\\fP [\\fIin2.fasta ... \\fP]");
  addDescription(parser, "Program for checksum of sequence reads. ");

  addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE,"fastafiles", true));

  setValidValues(parser, 0,"fa fa.gz fasta fasta.gz");

  addSection(parser, "Options");
  //add debug option:
  addOption(parser, seqan::ArgParseOption("d", "debug", "Debug mode. Prints full hex for each read to stdout"));
  addOption(parser, seqan::ArgParseOption("R", "no-readnames", "Do not use read names as part of checksum"));

  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res;
  }

  options.debug = seqan::isSet(parser, "debug");
  options.noReadNames = seqan::isSet(parser, "no-readnames");


  options.fastafiles = getArgumentValues(parser, 0);

  
  return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const **argv) {
  Fastainfo info; // Define structure variable
  seqan::ArgumentParser::ParseResult res = parseCommandLine(info, argc, argv); // Parse the command line.

  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res == seqan::ArgumentParser::PARSE_ERROR;
  }

  // Define:
  uint64_t sum = 0;
  unsigned count = 0;
  seqan::StringSet<seqan::CharString> idSub;
  seqan::CharString string2hash;

  seqan::CharString id;
  seqan::CharString seq;
  hash_t hex;

  // Open stream
  seqan::SeqFileIn seqFileIn;

  for (int i = 0; i < info.fastafiles.size(); i++) {
    const char* fasta = info.fastafiles[i].c_str();
    
    if (!open(seqFileIn, fasta)) {
      std::cerr << "ERROR: Could not open the file: " << fasta << " for reading.\n";
      return 1;
    }

    // Read record
    while (!seqan::atEnd(seqFileIn)) {
      try
      {
        readRecord(id, seq, seqFileIn);
      }
      catch (seqan::Exception const & e)
      {
        if (seqan::atEnd(seqFileIn)) {
          std::cerr << "WARNING: Could not continue reading " << fasta <<  " at line: " << count+1 << ".\n";
          return 1;
        }
        std::cerr << "ERROR: Could not read from " << fasta << "\n";
        return 1;
      }

      count +=1;

      // cut away after first space
      seqan::strSplit(idSub, id, seqan::EqualsChar<' '>(), false, 1);

      if (!info.noReadNames) {
        seqan::append(string2hash, idSub[0]);
        seqan::append(string2hash, "/1"); // to be consistent with BAM and FASTQ
      }
      seqan::append(string2hash, seq);

      // Get MD5 hash
      hex = str2md5(toCString(string2hash), length(string2hash));

      if (info.debug) {
        std::cout << string2hash << " " <<   std::hex << hex.p.low << "\n";
      } else {
        hexSum(hex, sum);
      }

      seqan::clear(string2hash);
      seqan::clear(idSub);
    }

  }

  if (!info.debug) {
    std::cout << std::hex << sum << "\t";
    std::cout << std::dec << count << "\n";
  }
    
  return 0;
}

