#include <iostream>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <stdint.h>

#include "bamhash_checksum_common.h"

static void usage(std::string name);
static void version(std::string name);


int main(int argc, char ** argv)
{
	
	int c;
	int debug = 0;
	while ((c = getopt (argc, argv, "vhd:")) != -1) {
		switch (c) {
		case 'h':
			// display help and exit normally
			usage(argv[0]);
			return 0;
		case 'v':
			// display version and exit normally
			version(argv[0]);
			return 0;	  
		case 'd':
			debug = 1;
			argv[1]=optarg;
			argv[2]=argv[3];
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	// Check arguments
	if (argc < 3 || argc > 4)
	{
		usage(argv[0]);
		return 0;
	}

	
	// Define:
	uint64_t sum = 0;
	unsigned count = 0;
	seqan::StringSet<seqan::CharString> idSub1;
	seqan::StringSet<seqan::CharString> idSub2;
	seqan::CharString string2hash1;
	seqan::CharString string2hash2;
	seqan::CharString id1;
	seqan::CharString id2;
	seqan::CharString seq1;
	seqan::CharString seq2;
	seqan::CharString qual1;
	seqan::CharString qual2;
	
	// Open GZStream
	seqan::Stream<seqan::GZFile> gzStream1;
	seqan::Stream<seqan::GZFile> gzStream2;
	
	if (!open(gzStream1, argv[1], "r")) {
		std::cerr << "ERROR: Could not open the file: " << argv[1] << " for reading.\n";
		return 1;
	}
	
	if (!open(gzStream2, argv[2], "r")) {
		std::cerr << "ERROR: Could not open the file: " << argv[2] << " for reading.\n";
		return 1;
	}
  
	//Setup RecordReader for reading FASTQ file from gzip-compressed file
	seqan::RecordReader<seqan::Stream<seqan::GZFile>, seqan::SinglePass<> > reader1(gzStream1);
	seqan::RecordReader<seqan::Stream<seqan::GZFile>, seqan::SinglePass<> > reader2(gzStream2);
  
	// Read record
	while (!atEnd(reader1) || !atEnd(reader2)) {
		
		if (readRecord(id1, seq1, qual1, reader1, seqan::Fastq()) != 0) {
	    if (atEnd(reader1)) {
				std::cerr << "WARNING: Could not continue reading " << argv[1] <<  " at line: " << count+1 << ". Check if files have the same number of reads.\n";
				return 1;
			}
	    std::cerr << "ERROR: Could not read from " << argv[1] << "\n";
	    return 1;
	  }
		
		if (readRecord(id2, seq2, qual2, reader2, seqan::Fastq()) != 0) {
	    if (atEnd(reader2)) {
				std::cerr << "WARNING: Could not continue reading " << argv[2] << " at line: " << count+1 << ". Check if files have the same number of reads.\n";
				return 1;
			}
	    std::cerr << "ERROR: Could not read from " << argv[2] << "\n";
	    return 1;
	  }
		
		count +=1;
		
		
		// If include id, then cut id on first whitespace
		if (seqan::endsWith(id1,"/1")) {
	    seqan::strSplit(idSub1, id1, '/', false, 1);
	  } else {
	    seqan::strSplit(idSub1, id1, ' ', false, 1);
	  }
	
		if (seqan::endsWith(id2,"/2")) {
	    seqan::strSplit(idSub2, id2, '/', false, 1);
	  } else {
	    seqan::strSplit(idSub2, id2, ' ', false, 1);
	  }

		// Check if names are in same order in both files
		if (!(idSub1[0] ==  idSub2[0]))
	  {
	    std::cerr << "WARNING: Id_names in line: " << count << " are not in the same order\n";
	    return 1;
	  }

		seqan::append(string2hash1, idSub1[0]);
		seqan::append(string2hash1,"/1");
		seqan::append(string2hash1, seq1);
		seqan::append(string2hash1, qual1);
		
		
		seqan::append(string2hash2, idSub2[0]);
		seqan::append(string2hash2,"/2");
		seqan::append(string2hash2, seq2);
		seqan::append(string2hash2, qual2);

		// Get MD5 hash
		hash_t hex1 = str2md5(toCString(string2hash1), length(string2hash1));
		hash_t hex2 = str2md5(toCString(string2hash2), length(string2hash2));

		if (debug !=0) {
	    std::cout << std::hex << hex1.p.low << "\n";
	    std::cout << std::hex << hex2.p.low << "\n";
	  } else {
	    hexSum(hex1, sum);
	    hexSum(hex2, sum);
	  }
		
		seqan::clear(string2hash1);
		seqan::clear(string2hash2);
		seqan::clear(idSub1);
		seqan::clear(idSub2);	

	}

	if (debug == 0) {
		std::cout << std::hex << sum << "\t";
		std::cout << std::dec << count << "\n";
	}
	
	return 0;
}

static void usage(std::string name)
{
	std::cout << "USAGE: " << name << "\t[-h|-v|-d]\n"
						<< "USAGE: " << name << "\t<in.R1.fastq.gz> <in.R2.fastq.gz>\n"
						<< "USAGE: " << name << "\t -d <in.R1.fastq.gz> <in.R2.fastq.gz>\n"
						<< "\t-h\tDisplay this help message\n"
						<< "\t-v\tDisplay version\n"
						<< "\t-d\tDebug mode. Prints full hex for each read to stdout\n";

}

static void version(std::string name)
{
	std::cout << "Version: " << BAMHASH_VERSION << std::endl;
}
