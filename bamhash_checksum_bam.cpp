#include <iostream>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <stdint.h>


#include "bamhash_checksum_common.h"


static void usage(std::string name);
static void version(std::string name);


/** only needed for seqan 1.4.1 and lower
inline bool
hasFlagSupplementary(seqan::BamAlignmentRecord const & record)
{
    return (record.flag & 0x0800) == 0x0800;
}
*/

int main(int argc, char ** argv)
{

    // Check arguments
    if (argc < 2 || argc > 3) {
      usage(argv[0]);
      return 0;
    }

    // Use getopt to read in optional arguments
    int c;
    int debug = 0;
    while ((c = getopt (argc, argv, "vhd:")) != -1) {
      switch (c)
			{
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
				break;
			default:
				usage(argv[0]);
				return 1;
			}
		}

    // Open BGZF Stream for reading.
    seqan::Stream<seqan::Bgzf> inStream;
    if (!open(inStream, argv[1], "r")) {
			std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
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
			std::cerr << "ERROR: Could not read header from BAM file " << argv[1] << "\n";
			return 1;
    }
    seqan::clear(header);
		
    // Define:
    seqan::BamAlignmentRecord record;
    uint64_t sum = 0;
    uint64_t count = 0;
    seqan::CharString string2hash;
    //char hexCstr[33];
		
    // Read record
    while (!atEnd(inStream)) {
      if (readRecord(record, context, inStream, seqan::Bam()) != 0) {
				std::cerr << "ERROR: Could not read record from BAM File " << argv[1] << "\n";
				return 1;
			}	
      count +=1;
      // Check if flag: reverse complement and change record accordingly
      if (hasFlagRC(record)) {
				reverseComplement(record.seq);
				reverse(record.qual);
			}
      // Check if flag: supplementary and exclude those
      if (!hasFlagSupplementary(record)) {
				// Construct one string from record
				seqan::append(string2hash, record.qName);
				if (hasFlagFirst(record)) {
					seqan::append(string2hash, "/1");
				} else if (hasFlagLast(record)) {
					seqan::append(string2hash, "/2");
				} else {
					std::cerr << "ERROR: A record in " << argv[1] << " has neither First or Last Flag.\n";
					return 1;
				}
				
				seqan::append(string2hash, record.seq);
				seqan::append(string2hash, record.qual);
				seqan::clear(record);
				
				// Get MD5 hash
				hash_t hex = str2md5(toCString(string2hash), length(string2hash));
				seqan::clear(string2hash);
				
				if (debug !=0) {
					std::cout << std::hex << hex.p.low << "\n";
				} else {
					hexSum(hex, sum); 
				}
			}
    }

    // print result
    if (debug == 0) {
			std::cout << std::hex << sum << "\t";
			std::cout << std::dec << count << "\n";
		}
		
    return 0;
}


static void usage(std::string name)
{
	std::cout << "USAGE: " << name << "\t[-h|-v|-d]\n"
						<< "USAGE: " << name << "\t<in.bam>\n"
						<< "USAFE: " << name << "\t -d <in.bam>\n"
						<< "\t-h\tDisplay this help message\n"
						<< "\t-v\tDisplay version\n"
						<< "\t-d\tDebug mode. Prints full hex for each read to stdout" << std::endl;;
	
}

static void version(std::string name)
{
	std::cout << "Version: " << BAMHASH_VERSION << std::endl;
}
