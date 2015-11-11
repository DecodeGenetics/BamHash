SHELL := /bin/bash

FASTQBIN=../bamhash_checksum_fastq
FASTABIN=../bamhash_checksum_fasta
BAMBIN=../bamhash_checksum_bam
WGSIM=wgsim
SEQTK=seqtk

REF=ecoliMG1655.fa

.PHONY: FORCE clean reallyclean


#r.namesorted.fastq.md5sum FORCE
r: r.sorted.bam.md5sum r.unsorted.bam.md5sum r.fastq.md5sum r.interleaved.md5sum r.namesorted.bam.md5sum r.single.sorted.bam.md5sum r.single.unsorted.bam.md5sum r.single.fastq.md5sum r.single.namesorted.bam.md5sum r.repeat.unsorted.bam.md5sum r.repeat.sorted.bam.md5sum r.repeat.mixed.bam.md5sum r.repeat.fastq.md5sum r.fasta.md5sum r.single.noqual.bam.md5sum r.single.noqual.fastq.md5sum FORCE
	for x in $^; do if [ $$x != "FORCE" ]; then echo $$x $$(cat $$x); fi; done

# retrieve ref from internet if not already available
${REF}:
	wget -nH https://github.com/ctb/edda/raw/master/doc/tutorials-2012/files/ecoliMG1655.fa.gz
	gunzip ${REF}.gz

${REF}.bwt: ${REF}
	bwa index ${REF}

r1.fastq: ${REF}.bwt
	${WGSIM} -N 100000 ${REF} r1.fastq r2.fastq


r2.fastq: r1.fastq
	#noop

r.interleaved.fastq: r1.fastq
	paste <(paste - - - - < r1.fastq) <(paste - - - - < r2.fastq)  | tr '\t' '\n'  > r.interleaved.fastq

r1.fasta: r1.fastq
	${SEQTK} seq -a r1.fastq | sed 's/\/[1-2]$$//' > r1.fasta

r.unsorted.bam: r1.fastq r2.fastq
	bwa mem ${REF} r1.fastq r2.fastq | samtools view -bSh - > r.unsorted.bam

r.single.unsorted.bam: r1.fastq
	bwa mem ${REF} r1.fastq | samtools view -bSh - > r.single.unsorted.bam


r.sorted.bam: r.unsorted.bam
  #samtools sort -o r.sorted.bam -T test -O bam r.unsorted.bam
	samtools sort  -f r.unsorted.bam r.sorted.bam

r.namesorted.bam: r.unsorted.bam
  #samtools sort -o r.namesorted.bam -T test -O bam -n r.unsorted.bam
	samtools sort -f -n r.unsorted.bam r.namesorted.bam



r.single.sorted.bam: r.single.unsorted.bam
	samtools sort  -f r.single.unsorted.bam r.single.sorted.bam

r.single.namesorted.bam: r.single.unsorted.bam
  #samtools sort -o r.namesorted.bam -T test -O bam -n r.unsorted.bam
	samtools sort -f -n r.single.unsorted.bam r.single.namesorted.bam


#r.namesorted.1.fastq: r.namesorted.bam
#	bam2fastq -p -o r.namesorted.fastq r.namesorted.bam



r.sorted.bam.md5sum: r.sorted.bam FORCE
	${BAMBIN} r.sorted.bam > r.sorted.bam.md5sum

r.namesorted.bam.md5sum: r.namesorted.bam FORCE
	${BAMBIN} r.namesorted.bam > r.namesorted.bam.md5sum

r.unsorted.bam.md5sum: r.unsorted.bam FORCE
	${BAMBIN} r.unsorted.bam > r.unsorted.bam.md5sum


r.single.sorted.bam.md5sum: r.single.sorted.bam FORCE
	${BAMBIN} -P r.single.sorted.bam > r.single.sorted.bam.md5sum

r.single.namesorted.bam.md5sum: r.single.namesorted.bam FORCE
	${BAMBIN} -P r.single.namesorted.bam > r.single.namesorted.bam.md5sum

r.single.unsorted.bam.md5sum: r.single.unsorted.bam FORCE
	${BAMBIN} r.single.unsorted.bam > r.single.unsorted.bam.md5sum

r.single.noqual.bam.md5sum: r.single.sorted.bam FORCE
	${BAMBIN} -Q r.single.sorted.bam > r.single.noqual.bam.md5sum

r.repeat.unsorted.bam.md5sum: r.unsorted.bam FORCE
	${BAMBIN} r.unsorted.bam r.unsorted.bam > r.repeat.unsorted.bam.md5sum

r.repeat.sorted.bam.md5sum: r.sorted.bam FORCE
	${BAMBIN} r.unsorted.bam r.unsorted.bam > r.repeat.sorted.bam.md5sum

r.repeat.mixed.bam.md5sum: r.unsorted.bam r.sorted.bam FORCE
	${BAMBIN} r.unsorted.bam r.sorted.bam > r.repeat.mixed.bam.md5sum


#r.namesorted.fastq.md5sum: r.namesorted.1.fastq FORCE
#	${FASTQBIN} r.namesorted.1.fastq r.namesorted.2.fastq > r.namesorted.fastq.md5sum

r.fastq.md5sum: r1.fastq r2.fastq FORCE
	${FASTQBIN} r1.fastq r2.fastq > r.fastq.md5sum

r.interleaved.md5sum: r.interleaved.fastq FORCE
	${FASTQBIN} --interleaved r.interleaved.fastq > r.interleaved.md5sum

r.single.fastq.md5sum: r1.fastq  FORCE
	${FASTQBIN} -P r1.fastq  > r.single.fastq.md5sum

r.single.noqual.fastq.md5sum: r1.fastq  FORCE
	${FASTQBIN} -P -Q r1.fastq  > r.single.noqual.fastq.md5sum


r.repeat.fastq.md5sum: r1.fastq r2.fastq FORCE
	${FASTQBIN} r1.fastq r2.fastq r1.fastq r2.fastq > r.repeat.fastq.md5sum


r.fasta.md5sum: r1.fasta FORCE
	${FASTABIN} r1.fasta > r.fasta.md5sum


FORCE:


clean:
	rm -f *.md5sum

reallyclean: clean
	rm -f *.bam *.fastq
