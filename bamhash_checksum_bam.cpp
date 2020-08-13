#include "htslib/hts.h"
#include <seqan/bam_io.h>
#include <iostream>
#include <seqan/stream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <seqan/arg_parse.h>
#include <seqan/hts_io.h>



#include "bamhash_checksum_common.h"

struct Baminfo {
  std::vector<std::string>  bamfiles;
  bool debug;
  bool noReadNames;
  bool noQuality;
  bool paired;
  seqan::CharString reference;

  Baminfo() : debug(false), noReadNames(false), noQuality(false), paired(true), reference("") {}

};

struct Counts {
  uint64_t sum;
  uint64_t count;

  Counts() : sum(0), count(0) {}

};

seqan::ArgumentParser::ParseResult
parseCommandLine(Baminfo& options, int argc, char const **argv) {
  // Setup ArgumentParser.
  seqan::ArgumentParser parser("bamhash_checksum_bam");
  //readlink("/proc/self/exe", options.bindir, sizeof(options.bindir)-1);

  setShortDescription(parser, "Checksum of a sam, bam or cram file"); //TODO change description
  setVersion(parser, BAMHASH_VERSION);
  setDate(parser, "Oct 2018");

  addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI<in.bam> <in2.bam> ...\\fP");
  addDescription(parser, "Program for checksum of sequence reads. ");

  addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE,"bamfile", "False", 1));
  setValidValues(parser, 0,"sam bam cram");

  addSection(parser, "Options");

  //add debug option:
  addOption(parser, seqan::ArgParseOption("d", "debug", "Debug mode. Prints full hex for each read to stdout"));
  addOption(parser, seqan::ArgParseOption("R", "no-readnames", "Do not use read names as part of checksum"));
  addOption(parser, seqan::ArgParseOption("Q", "no-quality", "Do not use read quality as part of checksum"));
  addOption(parser, seqan::ArgParseOption("P", "no-paired", "Cram files were not generated with paired-end reads"));
  addOption(parser, seqan::ArgParseOption("r", "reference-file", "Path to reference-file if reference not given in header",
                    seqan::ArgParseArgument::INPUT_FILE));

  setValidValues(parser, "reference-file", "fa");

  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res;
  }

  options.debug = isSet(parser, "debug");
  options.noReadNames = isSet(parser, "no-readnames");
  options.noQuality = isSet(parser, "no-quality");
  options.paired = !isSet(parser, "no-paired");
  getOptionValue(options.reference, parser, "reference-file");

  options.bamfiles = getArgumentValues(parser, 0);

  return seqan::ArgumentParser::PARSE_OK;
}


// -----------------------------------------------------------------------------
// FUNCTION getSampleIdAndLaneNames()
// -----------------------------------------------------------------------------

void getLaneNames(std::map<seqan::CharString, unsigned> & laneNames, std::string const & header)
{
  for (int i = 0; i < header.size(); /*empty on purpose*/)
  {
    auto hdr_find_it = std::find(header.begin() + i, header.end(), '\n');
    std::string line = header.substr(i, hdr_find_it - header.begin() - i);

    if (line.size() > 7 && line[0] == '@' && line[1] == 'R' && line[2] == 'G' && line[3] == '\t')
    {
      for (int j = 0; j < static_cast<int>(line.size()); /*empty on purpose*/)
      {
        auto line_find_it = std::find(line.begin() + j, line.end(), '\t');
        std::string field = line.substr(j, line_find_it - line.begin() - j);

        if (field.size() > 3 && field[0] == 'I' && field[1] == 'D' && field[2] == ':')
        {
          seqan::CharString read_group_name = field.substr(3);
          int read_group_index = laneNames.size();
          laneNames[read_group_name] = read_group_index;
        }

        j = std::distance(line.begin(), line_find_it) + 1;
      }
    }

    i = std::distance(header.begin(), hdr_find_it) + 1;
  }
}

// -----------------------------------------------------------------------------
// FUNCTION getLane()
// -----------------------------------------------------------------------------

int getLane(seqan::BamAlignmentRecord & record,
            seqan::BamTagsDict & tagsDict,
            std::map<seqan::CharString, unsigned> & laneNames)
{
  unsigned tagIdx = 0;

  if (!seqan::findTagKey(tagIdx, tagsDict, "RG"))
  {
    std::cerr << "ERROR: Found a read with a missing read group (RG) tag\n";
    return -1;
  }

  seqan::CharString read_group;

  if(!seqan::extractTagValue(read_group, tagsDict, tagIdx))
  {
    std::cerr << "ERROR: Failed to extract read group (RG) tag value\n";
    return -1;
  }

  return laneNames[read_group];
}



int main(int argc, char const **argv) {

  Baminfo info; // Define structure variable
  seqan::ArgumentParser::ParseResult res = parseCommandLine(info, argc, argv); // Parse the command line.

  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res == seqan::ArgumentParser::PARSE_ERROR;
  }

//Moving below to struct
//  uint64_t sum = 0;
//  uint64_t count = 0;
  bool pairedWarning = false;

  //adding new stuff
  std::map<seqan::CharString, unsigned> laneNames;
  // Initialize all counts for each lane.
  seqan::String<Counts> counts;

  for (int i = 0; i < info.bamfiles.size(); i++) {

    const char* bamfile = info.bamfiles[i].c_str();
    const char* reference = toCString(info.reference);

    // Open stream for reading
    seqan::HtsFile inStream(bamfile, "r", reference);

    // Initialize lane names (read groups).
    std::string header(inStream.hdr->text, inStream.hdr->l_text);
    getLaneNames(laneNames, header);
    unsigned lanecount = laneNames.size();
    resize(counts, lanecount);

    // Define:
    seqan::BamAlignmentRecord record;
    seqan::CharString string2hash;

    // Read record
    while (seqan::readRecord(record, inStream)){
      seqan::BamTagsDict tagsDict(record.tags);
      int l = getLane(record, tagsDict, laneNames);
      if (l == -1) return 1;

      // Check if flag: reverse complement and change record accordingly
      if (hasFlagRC(record)) {
        seqan::reverseComplement(record.seq);
        seqan::reverse(record.qual);
      }
      // Check if flag: supplementary and exclude those
      if (!hasFlagSupplementary(record) && !hasFlagSecondary(record)) {
        counts[l].count +=1;
        // Construct one string from record
        if (!info.noReadNames) {
          seqan::append(string2hash, record.qName);
          if(hasFlagLast(record)) {
            if (info.paired) {
              seqan::append(string2hash, "/2");
            } else {
              if (!pairedWarning) {
                std::cerr << "WARNING: seqread was run with --no-paired mode, but BAM file has reads marked as second pair" << std::endl;
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
          if (length(record.seq) == 1 && record.qual == " ")
            record.qual = "*";
          seqan::append(string2hash, record.qual);
        }
        seqan::clear(record);

        // Get MD5 hash
        hash_t hex = str2md5(toCString(string2hash), length(string2hash));

        if (info.debug) {
          std::cout << string2hash << " " << std::hex << hex.p.low << "\n";
        } else {
          hexSum(hex, counts[l].sum);
        }

        seqan::clear(string2hash);
      }
    }
  }

  if (!info.debug) {
    for (std::map<seqan::CharString, unsigned>::iterator it = laneNames.begin(); it != laneNames.end(); ++it) {
      std::cout << it->first << "\t";
      int lid = it->second;
      std::cout << std::hex << counts[lid].sum << "\t";
      std::cout << std::dec << counts[lid].count << "\n";
    }
  }

  return 0;
}
