# include SeqAn libraries, don't warn about MD5 deprecation
CXXFLAGS+=-I. -Wno-deprecated-declarations


# RELEASE build
CXXFLAGS+= -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1 
LDLIBS=-lz -lssl -lcrypto


TARGET = bamhash_checksum_bam bamhash_checksum_fastq bamhash_checksum_fasta
all: $(TARGET)

bamhash_checksum_bam: bamhash_checksum_common.o bamhash_checksum_bam.o
	 $(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bamhash_checksum_fastq: bamhash_checksum_common.o bamhash_checksum_fastq.o
	 $(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bamhash_checksum_fasta: bamhash_checksum_common.o bamhash_checksum_fasta.o
	 $(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)


clean:
	$(RM) *.o *~ $(TARGET)
