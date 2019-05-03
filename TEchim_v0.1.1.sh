#!/bin/bash

################################################################################
# TITLE: TEchim - detection of TE-gene chimera in RNA-seq data
# VERSION: 0.1.1 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 03/05/2019 (dd/mm/yyyy)
################################################################################


################################################################################
# Converting overlapping paired-end reads to in-silico pairs with 60nt on each #
# end, preserving the merged reads and preserving the strand orientation.      #
################################################################################

# set parameters
working_directory=/Users/koalcan/Documents/2019MAY_TEscoex_longRNA/TEST
INPUT_FASTQ1=test_1.fastq.gz
INPUT_FASTQ2=test_2.fastq.gz
SAMPLENAME=test
ncores=20
refgenome_dmel=/Users/koalcan/Documents/REF_2019FEB/dmel625_v03/
dmel_noTEs=/Users/koalcan/Documents/REF_2019FEB/dmel625_v04/dmel625_v04_noTEs.fa
dmel_onlyTEs=/Users/koalcan/Documents/REF_2019FEB/dmel625_v04/dmel625_v04_onlyTEs.fa

# change to wd
cd $working_directory

echo "starting"
date

# merge reads
flash $INPUT_FASTQ1 $INPUT_FASTQ2 \
	-z \
	-x 0.15 \
	-M 170 \
	-o $SAMPLENAME$"_out1" \
	-t $ncores

echo "done with merge"
date

################################################################################

# concatenate reads that were successfully merged with unmerged file _1.fastq
# the unmerged file _1.fastq will contain reads that are longer than 119nt,
# which can also be used to generate in-silico paired-end reads
cat $SAMPLENAME$"_out1.extendedFrags.fastq.gz" $SAMPLENAME$"_out1.notCombined_1.fastq.gz" > $SAMPLENAME$"_out1.combined.fastq.gz"

# for the following awk commands, the reads will be filtered. only reads
# containing at least 119nt are kept. this is to avoid overlap when taking 60nt
# sections from each end. 

# split up fasta into 4 separate files, one for each line. (a) header;
# (b) sequence; (c) +; (d) quality score.
# reads used here: merged (_m) _1.fastq.gz (_1) -> _m1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print a}}' <(gzip -dc $SAMPLENAME$"_out1.combined.fastq.gz") > a_m1.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print b}}' <(gzip -dc $SAMPLENAME$"_out1.combined.fastq.gz") > b_m1.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print c}}' <(gzip -dc $SAMPLENAME$"_out1.combined.fastq.gz") > c_m1.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print d}}' <(gzip -dc $SAMPLENAME$"_out1.combined.fastq.gz") > d_m1.1

# reads used here: _2.fastq.gz reads (2)
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print a}}' <(gzip -dc $SAMPLENAME$"_out1.notCombined_2.fastq.gz") > a_2.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print b}}' <(gzip -dc $SAMPLENAME$"_out1.notCombined_2.fastq.gz") > b_2.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print c}}' <(gzip -dc $SAMPLENAME$"_out1.notCombined_2.fastq.gz") > c_2.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print d}}' <(gzip -dc $SAMPLENAME$"_out1.notCombined_2.fastq.gz") > d_2.1

################################################################################

# add :A to readname of _m1 (This is to not confuse with otherwise
# duplicate read names of _2)
sed 's/ /:A /g' a_m1.1 > a_m1.2
# crop sequences (60nt) from _m1 which will be used for NEW_1.fasta
sed -E 's/(.{60}).*/\1/' b_m1.1 > b_m1.2
# crop quality scores (60) from _m1 which will be used for NEW_1.fasta
sed -E 's/(.{60}).*/\1/' d_m1.1 > d_m1.2

# crop sequences (60nt) and generate reverse-complement from _m1 which will be
# used for NEW_2.fasta
rev b_m1.1 | sed -E 's/(.{60}).*/\1/' | grep '^[ATCGN]' /dev/stdin | tr ATCGN TAGCN > b_m1.3
# crop quality scores (60) and invert from _m1 which will be used for
# NEW_2.fasta
rev d_m1.1 | sed -E 's/(.{60}).*/\1/' > d_m1.3

# combine in-silico components for NEW_1.fasta file
paste -d'\n' a_m1.2 b_m1.2 c_m1.1 d_m1.2 > $SAMPLENAME$"_out3_1.fasta"
# combine in-silico components for NEW_2.fasta file
paste -d'\n' a_m1.2 b_m1.3 c_m1.1 d_m1.3 > $SAMPLENAME$"_out3_2.fasta"

################################################################################

# add :B to readname of _2 (This is to not confuse with otherwise
# duplicate read names of _m1)
sed 's/ /:B /g' a_2.1 > a_2.2
# crop sequences (60nt) from _2 which will be used for NEW_2.fasta
sed -E 's/(.{60}).*/\1/' b_2.1 > b_2.2
# crop quality scores (60) from _2 which will be used for NEW_2.fasta
sed -E 's/(.{60}).*/\1/' d_2.1 > d_2.2

# crop sequences (60nt) and generate reverse-complement from _2 which will be
# used for NEW_1.fasta
rev b_2.1 | sed -E 's/(.{60}).*/\1/' | grep '^[ATCGN]' /dev/stdin | tr ATCGN TAGCN > b_2.3
# crop quality scores (60) and invert from _2 which will be used for
# NEW_1.fasta
rev d_2.1 | sed -E 's/(.{60}).*/\1/' > d_2.3

# combine in-silico components and add to NEW_1.fasta file
paste -d'\n' a_2.2 b_2.3 c_2.1 d_2.3 >> $SAMPLENAME$"_out3_1.fasta"
# combine in-silico components and add to NEW_2.fasta file
paste -d'\n' a_2.2 b_2.2 c_2.1 d_2.2 >> $SAMPLENAME$"_out3_2.fasta"

################################################################################

# Create a lookup table for long reads. The long read is the ACTUAL nucleotide
# sequence of the mRNA molecule (i.e. the "READ_2 from the stranded RNA library
# kit)

# combine _m1 and _2 read names. use only the read name (not the tab-separated
# additional flag)
cat a_m1.2 a_2.2 | sed 's/ .*//' > a_m12.1

# create reverse-complement of the _m1 sequences. this means that the long reads
# will represent the actual nucleotide sequence of the mRNA strand
rev b_m1.1 | grep '^[ATCGN]' /dev/stdin | tr ATCGN TAGCN > b_m1.4

# combine _m1 and _2 sequences
cat b_m1.4 b_2.1 > b_m12.1

# create look-up table and remove "@" at the beginning of the readname
paste -d'\t' a_m12.1 b_m12.1 | sed 's/^@\(.*\)/\1/' > $SAMPLENAME$"_LOOKUP.tsv"

rm a_*
rm b_*
rm c_*
rm d_*

echo "done with cropping"
date

################################################################################

# run STAR in chimera-mode
STAR --runThreadN $ncores \
	--genomeDir $refgenome_dmel \
	--readFilesIn $SAMPLENAME$"_out3_1.fasta" $SAMPLENAME$"_out3_2.fasta" \
	--chimSegmentMin 20 \
	--chimOutType WithinBAM \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix $SAMPLENAME$"_out4_" \
	&&
	
# extract only hits that cross TE-GENE breakpoints. the awk commands remove
# (1) TE-TE reads and (2)&(3) reads that span the TE and it's LTR
samtools view $SAMPLENAME$"_out4_"Aligned.sortedByCoord.out.bam | grep TEchr_ | awk '($7 != "=") && ($3 !~ $7) && ($7 !~ $3)' > $SAMPLENAME$"_out5_TExGENES.sam"

echo "done with mapping"
date

# extract readnames and remove duplicates
awk '{print $1}' $SAMPLENAME$"_out5_TExGENES.sam" | sort | uniq > $SAMPLENAME$"_out6_TExGENES_readnames.txt"

# combine readnames with long sequences stored in the lookup file
join -1 1 -2 1 <(sort $SAMPLENAME$"_out6_TExGENES_readnames.txt") <(sort $SAMPLENAME$"_LOOKUP.tsv") > $SAMPLENAME$"_out7_TExGENES_longreads.tsv"

################################################################################

# the following while loop will add data to file using ">>". just as safety
# precaution, this line makes sure no data with such a name exists
rm -f $SAMPLENAME$"_out8_TExGENES_blastedreads_plusnohit.tsv"

# loop through long sequences. two blast searches are performed, one on a
# reference genome without TE sequences and another one on all TE sequences.
# only the top hit is reported.
while IFS=$' ' read -r readname sequence
do
		echo "${readname}" > var1
		echo "${sequence}" > var2
		echo "${sequence}" | blastn -db $dmel_noTEs -outfmt "6 qstart qend sseqid sstart send sstrand" -num_alignments 1 -num_threads $ncores | head -n 1 > var3		
		if ! [ -s var3 ]; then echo "no-hit" > var3; fi
		echo "${sequence}" | blastn -db $dmel_onlyTEs -outfmt "6 qstart qend sseqid sstart send slen sstrand" -num_alignments 1 -num_threads $ncores | head -n 1 > var4
		if ! [ -s var4 ]; then echo "no-hit" > var4; fi
		paste var1 var2 var3 var4 >> $SAMPLENAME$"_out8_TExGENES_blastedreads_plusnohit.tsv"
		rm var*
done < $SAMPLENAME$"_out7_TExGENES_longreads.tsv"

# remove reads that did not give BLAST result
grep -v no-hit $SAMPLENAME$"_out8_TExGENES_blastedreads_plusnohit.tsv" > $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv"

echo "done with BLAST"
date

################################################################################

# determine whether the section that maps to the genome is the 5' or the 3' end
# of the mRNA section:
# 5'-######|GENE|TE|######-3' => PART5
# 5'-######|TE|GENE|######-3' => PART3
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; if (a < b) {print "PART5"} else {print "PART3"}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.genepart

# determine the precise breakpoint on the chromosome. this depends on  whether
# the chromosomal part is PART5 or PART3
awk 'BEGIN {OFS = "\t"} {a = $5; {print a}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.chr
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d-1} else {print c-1}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.start
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d} else {print c}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.end
awk 'BEGIN {OFS = "\t"} {a = $1 ; {print a}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.readname
awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "."}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.score
awk 'BEGIN {OFS = "\t"} {a = $8 ; if (a == "plus") {print "+"} else {print "-"}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.strand

paste -d'\t' tmpfile.chr tmpfile.breakpoint.chr.start tmpfile.breakpoint.chr.end tmpfile.readname tmpfile.score tmpfile.breakpoint.chr.strand > $SAMPLENAME$"out10_breakpoints.bed"

# determine the precise breakpoint on the TE. this depends on  whether the TE
# part is PART5 or PART3
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $13 ; d = $12 ; if (a < b) {print d} else {print c}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.TE

# determine the overlap between the two mapped sections of the long read
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $4 ; d = $10 ; if (a < b) {print b-c-1} else {print a-d-1}}' < $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" > tmpfile.uncertainty

paste -d'\t' $SAMPLENAME$"_out9_TExGENES_blastedreads.tsv" tmpfile.breakpoint.chr.end tmpfile.breakpoint.TE tmpfile.uncertainty > $SAMPLENAME$"out11_combined_results.tsv"

echo "done!"
date

