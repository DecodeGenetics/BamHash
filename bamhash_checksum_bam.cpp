#include <iostream>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <seqan/arg_parse.h>


#include "bamhash_checksum_common.h"


/** only needed for seqan 1.4.1 and lower
inline bool
hasFlagSupplementary(seqan::BamAlignmentRecord const & record)
{
    return (record.flag & 0x0800) == 0x0800;
}
*/
struct Baminfo {
  std::vector<std::string>  bamfiles;
  bool debug;
  bool noReadNames;
  bool noQuality;
  bool paired;

  Baminfo() : debug(false), noReadNames(false), noQuality(false), paired(true) {}

};

seqan::ArgumentParser::ParseResult
parseCommandLine(Baminfo& options, int argc, char const **argv) {
  // Setup ArgumentParser.
  seqan::ArgumentParser parser("bamhash_checksum_bam");
  //readlink("/proc/self/exe", options.bindir, sizeof(options.bindir)-1);

  setShortDescription(parser, "Checksum of a bam file"); //TODO change description
  setVersion(parser, BAMHASH_VERSION);
  setDate(parser, "May 2015");

  addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI<in.bam> <in2.bam> ...\\fP");
  addDescription(parser, "Program for checksum of sequence reads. ");

  addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE,"bamfile", "False", 1));

  setValidValues(parser, 0,"bam sam");

  addSection(parser, "Options");

  //add debug option:
  addOption(parser, seqan::ArgParseOption("d", "debug", "Debug mode. Prints full hex for each read to stdout"));
  addOption(parser, seqan::ArgParseOption("R", "no-readnames", "Do not use read names as part of checksum"));
  addOption(parser, seqan::ArgParseOption("Q", "no-quality", "Do not use read quality as part of checksum"));
  addOption(parser, seqan::ArgParseOption("P", "no-paired", "Bam files were not generated with paired-end reads"));

  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res;
  }

  options.debug = isSet(parser, "debug");
  options.noReadNames = isSet(parser, "no-readnames");
  options.noQuality = isSet(parser, "no-quality");
  options.paired = !isSet(parser, "no-paired");

  options.bamfiles = getArgumentValues(parser, 0);

  return seqan::ArgumentParser::PARSE_OK;
}


int main(int argc, char const **argv) {

  Baminfo info; // Define structure variable
  seqan::ArgumentParser::ParseResult res = parseCommandLine(info, argc, argv); // Parse the command line.

  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res == seqan::ArgumentParser::PARSE_ERROR;
  }

  uint64_t sum = 0;
  uint64_t count = 0;
  bool pairedWarning = false;

  

  for (int i = 0; i < info.bamfiles.size(); i++) {
    // Open BGZF Stream for reading.
    seqan::Stream<seqan::Bgzf> inStream;

    const char* bamfile = info.bamfiles[i].c_str();
    
    if (!open(inStream, bamfile, "r")) {
      std::cerr << "ERROR: Could not open " << bamfile << " for reading.\n";
      return 1;
    }

    // Setup name store, cache, and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);

    // Read header.
    seqan::BamHeader header;
    if (readRecord(header, context, inStream, seqan::Bam()) != 0) {
      std::cerr << "ERROR: Could not read header from BAM file " << bamfile << "\n";
      return 1;
    }
    seqan::clear(header);

    // Define:
    seqan::BamAlignmentRecord record;
    seqan::CharString string2hash;
    //char hexCstr[33];

    // Read record
    while (!atEnd(inStream)) {
      if (readRecord(record, context, inStream, seqan::Bam()) != 0) {
        std::cerr << "ERROR: Could not read record from BAM File " << bamfile << "\n";
        return 1;
      }
      // Check if flag: reverse complement and change record accordingly
      if (hasFlagRC(record)) {
        reverseComplement(record.seq);
        reverse(record.qual);
      }
      // Check if flag: supplementary and exclude those
      if (!hasFlagSupplementary(record) && !hasFlagSecondary(record)) {
        count +=1;
        // Construct one string from record
        if (!info.noReadNames) {
          seqan::append(string2hash, record.qName);
          if(hasFlagLast(record)) {
            if (info.paired) {
              seqan::append(string2hash, "/2");
            } else {
              if (!pairedWarning) {
                std::cerr << "WARNING: BamHash was run with --no-paired mode, but BAM file has reads marked as second pair" << std::endl;
                pairedWarning = true;
              }
              seqan::append(string2hash, "/1");
            }
          } else {
            seqan::append(string2hash, "/1");
          }
        }

        seqan::append(string2hash, record.seq);
        if (!info.noQuality) {
          seqan::append(string2hash, record.qual);
        }
        seqan::clear(record);

        // Get MD5 hash
        hash_t hex = str2md5(toCString(string2hash), length(string2hash));


        if (info.debug) {
          std::cout << string2hash << " " << std::hex << hex.p.low << "\n";
        } else {
          hexSum(hex, sum);
        }

        seqan::clear(string2hash);
      }
    }

    // print result
    
  }

  if (!info.debug) {
    std::cout << std::hex << sum << "\t";
    std::cout << std::dec << count << "\n";
  }
  
  return 0;
}

