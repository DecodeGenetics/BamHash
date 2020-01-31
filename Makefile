# include SeqAn libraries, don't warn about MD5 deprecation
CXXFLAGS+=-I. -Wno-deprecated-declarations

#include htslib by setting -I<path_to_htslib/include> and link to the htslib library
HTSDIR=<path_to_htsdir>
CXXFLAGS+=-I$(HTSDIR)/include

# RELEASE build
CXX=g++ -std=c++11 -pthread
CXXFLAGS+= -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1
LDFLAGS=-L$(HTSDIR)/lib -lz -lssl -lcrypto -Wl,-rpath,$(HTSDIR)/lib -lhts

TARGET = bamhash_checksum_bam bamhash_checksum_fastq bamhash_checksum_fasta
all: $(TARGET)

bamhash_checksum_bam: bamhash_checksum_common.o bamhash_checksum_bam.o
	 $(CXX) $(LDFLAGS) -o $@ $^

bamhash_checksum_fastq: bamhash_checksum_common.o bamhash_checksum_fastq.o
	 $(CXX) $(LDFLAGS) -o $@ $^

bamhash_checksum_fasta: bamhash_checksum_common.o bamhash_checksum_fasta.o
	 $(CXX) $(LDFLAGS) -o $@ $^

clean:
	$(RM) *.o *~ $(TARGET)
