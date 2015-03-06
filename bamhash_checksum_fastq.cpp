#include <iostream>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <stdint.h>
#include <seqan/arg_parse.h>

#include "bamhash_checksum_common.h"

struct Fastqinfo
{
    seqan::CharString fastq1;
    seqan::CharString fastq2;
    bool debug;

    Fastqinfo() :
        debug(false)
    {}

};

seqan::ArgumentParser::ParseResult
parseCommandLine(Fastqinfo & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bamhash_checksum_fastq");
    //readlink("/proc/self/exe", options.bindir, sizeof(options.bindir)-1);

    setShortDescription(parser, "Checksum of a set of fastq files");
    setVersion(parser, BAMHASH_VERSION);
    setDate(parser, "Feb 2015");

    //addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI<in1.fastq.gz>\\fP \\fI<in2.fastq.gz>\\fP");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI<in.fastq.gz>\\fP");
    addDescription(parser, "Program for checksum of sequence reads. ");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE,"fastqfile_1"));
    //addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE,"fastqfile_2"));

    setValidValues(parser, 0,"fastq fastq.gz");
    //setValidValues(parser, 1,"fastq fastq.gz");

    addSection(parser, "Options");
    //add debug option:
    addOption(parser, seqan::ArgParseOption("d", "debug", "Debug mode. Prints full hex for each read to stdout"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.debug = isSet(parser, "debug");
    getArgumentValue(options.fastq1, parser, 0);
    //getArgumentValue(options.fastq2, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;
}


int main(int argc, char const ** argv)
{
	Fastqinfo info; // Define structure variable
    seqan::ArgumentParser::ParseResult res = parseCommandLine(info, argc, argv); // Parse the command line.

    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }

	// Define:
	uint64_t sum = 0;
	unsigned count = 0;
	seqan::StringSet<seqan::CharString> idSub1;
//	seqan::StringSet<seqan::CharString> idSub2;
	seqan::CharString string2hash1;
//	seqan::CharString string2hash2;
	seqan::CharString id1;
//	seqan::CharString id2;
	seqan::CharString seq1;
//	seqan::CharString seq2;
	seqan::CharString qual1;
//	seqan::CharString qual2;
	
	// Open GZStream
	seqan::Stream<seqan::GZFile> gzStream1;
//	seqan::Stream<seqan::GZFile> gzStream2;
	
	if (!open(gzStream1, toCString(info.fastq1), "r")) {
		std::cerr << "ERROR: Could not open the file: " << info.fastq1 << " for reading.\n";
		return 1;
	}
	
	//if (!open(gzStream2, toCString(info.fastq2), "r")) {
	//	std::cerr << "ERROR: Could not open the file: " << info.fastq2 << " for reading.\n";
	//	return 1;
	//}
  
	//Setup RecordReader for reading FASTQ file from gzip-compressed file
	seqan::RecordReader<seqan::Stream<seqan::GZFile>, seqan::SinglePass<> > reader1(gzStream1);
	//seqan::RecordReader<seqan::Stream<seqan::GZFile>, seqan::SinglePass<> > reader2(gzStream2);
  
	// Read record
	//while (!atEnd(reader1) || !atEnd(reader2)) {
	while (!atEnd(reader1)) {
	    if (readRecord(id1, seq1, qual1, reader1, seqan::Fastq()) != 0) {
	        if (atEnd(reader1)) {
    	            std::cerr << "WARNING: Could not continue reading " << info.fastq1 <<  " at line: " << count+1 << ". Check if files have the same number of reads.\n";
                    return 1;
                }
	        std::cerr << "ERROR: Could not read from " << info.fastq1 << "\n";
	        return 1;
	    }
		
/*
	    if (readRecord(id2, seq2, qual2, reader2, seqan::Fastq()) != 0) {
	    if (atEnd(reader2)) {
            std::cerr << "WARNING: Could not continue reading " << info.fastq2 << " at line: " << count+1 << ". Check if files have the same number of reads.\n";
            return 1;
        }
	    std::cerr << "ERROR: Could not read from " << info.fastq2 << "\n";
	    return 1;
	  }
*/
		
		count +=1;
		
		
		// If include id, then cut id on first whitespace
		if (seqan::endsWith(id1,"/1")) {
	        seqan::strSplit(idSub1, id1, '/', false, 1);
	    } else {
	        seqan::strSplit(idSub1, id1, ' ', false, 1);
	    }
	
/*
		if (seqan::endsWith(id2,"/2")) {
	        seqan::strSplit(idSub2, id2, '/', false, 1);
	    } else {
	        seqan::strSplit(idSub2, id2, ' ', false, 1);
	    }
*/

		// Check if names are in same order in both files
/*
		if (!(idSub1[0] ==  idSub2[0]))
	    {
	        std::cerr << "WARNING: Id_names in line: " << count << " are not in the same order\n";
	        return 1;
	    }
*/

		seqan::append(string2hash1, idSub1[0]);
		seqan::append(string2hash1,"/1");
		seqan::append(string2hash1, seq1);
		seqan::append(string2hash1, qual1);
		
		
/*
		seqan::append(string2hash2, idSub2[0]);
		seqan::append(string2hash2,"/2");
		seqan::append(string2hash2, seq2);
		seqan::append(string2hash2, qual2);
*/

		// Get MD5 hash
		hash_t hex1 = str2md5(toCString(string2hash1), length(string2hash1));
//		hash_t hex2 = str2md5(toCString(string2hash2), length(string2hash2));

		if (info.debug) {
	        std::cout << std::hex << hex1.p.low << "\n";
//	        std::cout << std::hex << hex2.p.low << "\n";
	    } else {
	        hexSum(hex1, sum);
//	        hexSum(hex2, sum);
	    }
		
		seqan::clear(string2hash1);
//		seqan::clear(string2hash2);
		seqan::clear(idSub1);
//		seqan::clear(idSub2);	

	}

	if (!info.debug) {
		std::cout << std::hex << sum << "\t";
		std::cout << std::dec << count << "\n";
	}
	
	return 0;
}

