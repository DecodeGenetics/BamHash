# BamHash


Hash BAM and FASTQ files to verify data integrity

For each pair of reads in a BAM or FASTQ file we compute a hash value
composed of the readname, whether it is first or last in pair, sequence and quality value.
All the hash values are summed up so the result is independent of the ordering within the files.
The result can be compared to verify that the pair of FASTQ files contain the same read 
information as the aligned BAM file.

## Manuscript

Arna Óskarsdóttir, Gísli Másson and Páll Melsted (2015) 
BamHash: a checksum program for verifying the integrity of sequence data. 
[Bioinformatics, btv539](http://bioinformatics.oxfordjournals.org/content/early/2015/10/01/bioinformatics.btv539.abstract).

A preprint is available on [bioRxiv](http://biorxiv.org/content/early/2015/03/03/015867).

## Usage

The program has three executables which are used for different filetypes. Running them with `--help` displays detailed help messages.

### Common options

All programs work with sets of reads. The reads are made up of a read name, sequence and quality information. All of these components go into the hash, but the read name or quality information can be ignored if necessary. This would be the case if a pipeline mangled the names, quantizised the quality or after realigning quality scores.

The default mode is to assume paired end reads. If you have single end reads you can supply the `--no-paired` option.

A debug option `-d` prints the information and hash value of each read individually, this can be helpful if BamHash is not cooperating with your pipeline.

Both multiline FASTA and FASTQ are supported and gzipped input for FASTA and FASTQ.

### BAM

~~~
bamhash_checksum_bam [OPTIONS] <in.bam> <in2.bam> ...
~~~

processes a number of BAM files. BAM files are assumed to contain paired end reads. If you run with `--no-paired` it treats all reads as single end and displays a warning if any read is marked as "second in pair" in the BAM file.

### FASTQ

~~~
bamhash_checksum_fastq [OPTIONS] <in1.fastq.gz> [in2.fastq.gz ... ]
~~~

processes a number of FASTQ files. FASTQ files are assumed to contain paired end reads, such that the first two files contain the first pair of reads, etc. If any of the read names in the two pairs don't match the program exits with failure.

### FASTA

~~~
bamhash_checksum_fasta [OPTIONS] <in1.fasta> [in2.fasta ... ]
~~~

processes a number of FASTA files. All FASTA files are assumed to be single end reads with no quality information. To compare to a BAM file, run `bamhash_checksum_bam --no-paired --no-quality`

## Compiling

The only external dependency is on OpenSSL for the MD5 implementation.
