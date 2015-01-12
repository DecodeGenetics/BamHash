BamHash
=======

Hash BAM and FASTQ files to verify data integrity

For each pair of reads in a BAM or FASTQ file we compute a hash value
composed of the readname, whether it is first or last in pair, sequence and quality value.
All the hash values are summed up so the result is independent of the ordering within the files.
The result can be compared to verify that the pair of FASTQ files contain the same read 
information as the aligned BAM file.

Manuscript
==========

In preperation.

Compiling
=========

The only external dependency is on OpenSSL for the MD5 implementation.
