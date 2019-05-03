#!/bin/bash

################################################################################
# TITLE: TEchim - detection of TE-gene chimera in RNA-seq data
# VERSION: 0.1.1 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 03/05/2019 (dd/mm/yyyy)
# DESCRIPTION: This tool converts overlapping paired-end reads to in-silico
# pairs of 60nt length, which are then screened for pairs where one mate maps
# in the genome and the other on a transposon. The contigs of these read pairs
# are then used to determine the precise breakpoint in the genome and the TE.
################################################################################

################################################################################
# REQUIREMENTS:
# - flash (https://ccb.jhu.edu/software/FLASH/)
# - STAR (https://github.com/alexdobin/STAR)
# - blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
################################################################################

# set parameters
wd=~							# working directory
FASTQ1=READS_1.fastq.gz			# input FASTQ file 1 (FULL PATH)
FASTQ2=READS_2.fastq.gz			# input FASTQ file 2 (FULL PATH)
SNa=EXP							# sample name	
SNo=1							# sample number
LNo=1							# sequencing lane number
nc=20							# number of cores
refgenome_dmel=/Users/koalcan/Documents/REF_2019FEB/dmel625_v03/
dmel_noTEs=/Users/koalcan/Documents/REF_2019FEB/dmel625_v04/dmel625_v04_noTEs.fa
dmel_onlyTEs=/Users/koalcan/Documents/REF_2019FEB/dmel625_v04/dmel625_v04_onlyTEs.fa

# change to wd
cd $wd

mkdir $SNa"_S"$SNo"_L"$LNo
cd $SNa"_S"$SNo"_L"$LNo

echo "TEchim v 0.1.1 run" > $SNa"_S"$SNo"_L"$LNo".log"
echo "input fastq_1:" "$FASTQ1" >> $SNa"_S"$SNo"_L"$LNo".log"
echo "input fastq_2:" "$FASTQ2" >> $SNa"_S"$SNo"_L"$LNo".log"
echo "sample name:" "$SNa" >> $SNa"_S"$SNo"_L"$LNo".log"
echo "sample number:" "$SNo" >> $SNa"_S"$SNo"_L"$LNo".log"
echo "sequencing lane number:" "$SNo" >> $SNa"_S"$SNo"_L"$LNo".log"
echo "--------------------------------" >> $SNa"_S"$SNo"_L"$LNo".log"
echo " --> starting at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo".log"

# merge reads
flash $FASTQ1 $FASTQ2 \
	-z \
	-x 0.15 \
	-M 170 \
	-o $SNa"_S"$SNo"_L"$LNo$"_out1" \
	-t $nc

echo " --> done with merging at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo".log"

################################################################################

# concatenate reads that were successfully merged with unmerged file _1.fastq
# the unmerged file _1.fastq will contain reads that are longer than 119nt,
# which can also be used to generate in-silico paired-end reads
cat $SNa"_S"$SNo"_L"$LNo$"_out1.extendedFrags.fastq.gz" $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_1.fastq.gz" > $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz"

# for the following awk commands, the reads will be filtered. only reads
# containing at least 119nt are kept. this is to avoid overlap when taking 60nt
# sections from each end. 

# split up fasta into 4 separate files, one for each line. (a) header;
# (b) sequence; (c) +; (d) quality score.
# reads used here: merged (_m) _1.fastq.gz (_1) -> _m1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print a}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz") > a_m1.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print b}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz") > b_m1.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print c}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz") > c_m1.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print d}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz") > d_m1.1

# reads used here: _2.fastq.gz reads (2)
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print a}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_2.fastq.gz") > a_2.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print b}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_2.fastq.gz") > b_2.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print c}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_2.fastq.gz") > c_2.1
awk 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= 120) {print d}}' <(gzip -dc $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_2.fastq.gz") > d_2.1 &&

rm $SNa"_S"$SNo"_L"$LNo$"_out1"*

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
paste -d'\n' a_m1.2 b_m1.2 c_m1.1 d_m1.2 > $SNa"_S"$SNo"_L"$LNo$"_out3_1.fasta"
# combine in-silico components for NEW_2.fasta file
paste -d'\n' a_m1.2 b_m1.3 c_m1.1 d_m1.3 > $SNa"_S"$SNo"_L"$LNo$"_out3_2.fasta"

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
paste -d'\n' a_2.2 b_2.3 c_2.1 d_2.3 >> $SNa"_S"$SNo"_L"$LNo$"_out3_1.fasta"
# combine in-silico components and add to NEW_2.fasta file
paste -d'\n' a_2.2 b_2.2 c_2.1 d_2.2 >> $SNa"_S"$SNo"_L"$LNo$"_out3_2.fasta"

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
paste -d'\t' a_m12.1 b_m12.1 | sed 's/^@\(.*\)/\1/' > $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.tsv" &&

rm a_*
rm b_*
rm c_*
rm d_*

echo " --> done with cropping at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo".log"

################################################################################

# run STAR in chimera-mode
STAR --runThreadN $nc \
	--genomeDir $refgenome_dmel \
	--readFilesIn $SNa"_S"$SNo"_L"$LNo$"_out3_1.fasta" $SNa"_S"$SNo"_L"$LNo$"_out3_2.fasta" \
	--chimSegmentMin 20 \
	--chimOutType WithinBAM \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix $SNa"_S"$SNo"_L"$LNo$"_out4_"
	
# extract only hits that cross TE-GENE breakpoints. the awk commands remove
# (1) TE-TE reads and (2)&(3) reads that span the TE and it's LTR
samtools view $SNa"_S"$SNo"_L"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | grep TEchr_ | awk '($7 != "=") && ($3 !~ $7) && ($7 !~ $3)' > $SNa"_S"$SNo"_L"$LNo$"_out5_TExGENES.sam" &&

rm $SNa"_S"$SNo"_L"$LNo$"_out3"*

mkdir $SNa"_S"$SNo"_L"$LNo"_STAR"
mv $SNa"_S"$SNo"_L"$LNo$"_out4"* $SNa"_S"$SNo"_L"$LNo"_STAR"/.

echo " --> done with mapping at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo".log"

# extract readnames and remove duplicates
awk '{print $1}' $SNa"_S"$SNo"_L"$LNo$"_out5_TExGENES.sam" | sort | uniq > $SNa"_S"$SNo"_L"$LNo$"_out6_TExGENES_readnames.txt"

# combine readnames with long sequences stored in the lookup file
join -1 1 -2 1 <(sort $SNa"_S"$SNo"_L"$LNo$"_out6_TExGENES_readnames.txt") <(sort $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.tsv") > $SNa"_S"$SNo"_L"$LNo$"_out7_TExGENES_longreads.tsv" &&

rm $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.tsv"
rm $SNa"_S"$SNo"_L"$LNo$"_out6_TExGENES_readnames.txt"

################################################################################

# the following while loop will add data to file using ">>". just as safety
# precaution, this line makes sure no data with such a name exists
rm -f $SNa"_S"$SNo"_L"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"

# loop through long sequences. two blast searches are performed, one on a
# reference genome without TE sequences and another one on all TE sequences.
# only the top hit is reported.
while IFS=$' ' read -r readname sequence
do
		echo "${readname}" > var1
		echo "${sequence}" > var2
		echo "${sequence}" | blastn -db $dmel_noTEs -outfmt "6 qstart qend sseqid sstart send sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var3		
		if ! [ -s var3 ]; then echo "no-hit" > var3; fi
		echo "${sequence}" | blastn -db $dmel_onlyTEs -outfmt "6 qstart qend sseqid sstart send slen sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var4
		if ! [ -s var4 ]; then echo "no-hit" > var4; fi
		paste var1 var2 var3 var4 >> $SNa"_S"$SNo"_L"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
		rm var*
done < $SNa"_S"$SNo"_L"$LNo$"_out7_TExGENES_longreads.tsv"

# remove reads that did not give BLAST result
grep -v no-hit $SNa"_S"$SNo"_L"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv" > $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" &&

rm $SNa"_S"$SNo"_L"$LNo$"_out7_TExGENES_longreads.tsv"
rm $SNa"_S"$SNo"_L"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"

echo " --> done with BLAST at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo".log"

################################################################################

# determine whether the section that maps to the genome is the 5' or the 3' end
# of the mRNA section:
# 5'-######|GENE|TE|######-3' => PART5
# 5'-######|TE|GENE|######-3' => PART3
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; if (a < b) {print "PART5"} else {print "PART3"}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.genepart

# determine the precise breakpoint on the chromosome. this depends on  whether
# the chromosomal part is PART5 or PART3
awk 'BEGIN {OFS = "\t"} {a = $5; {print a}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.chr
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d-1} else {print c-1}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.start
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d} else {print c}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.end

awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
	a = $11
	gsub(/TEchr_/,"",a)
	b = $15
	c = $3
	d = $9
	if (c < d) {print a"|"b"|GENE-TE|S"s"|L"l} else {print a"|"b"|TE-GENE|S"s"|L"l}
	}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.readname


awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "."}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.score
awk 'BEGIN {OFS = "\t"} {a = $8 ; if (a == "plus") {print "+"} else {print "-"}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.strand

paste -d'\t' tmpfile.chr tmpfile.breakpoint.chr.start tmpfile.breakpoint.chr.end tmpfile.readname tmpfile.score tmpfile.breakpoint.chr.strand > $SNa"_S"$SNo"_L"$LNo$"_out10_breakpoints.bed"

# determine the precise breakpoint on the TE. this depends on  whether the TE
# part is PART5 or PART3
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $13 ; d = $12 ; if (a < b) {print d} else {print c}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.TE

# determine the overlap between the two mapped sections of the long read
awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $4 ; d = $10 ; if (a < b) {print b-c-1} else {print a-d-1}}' < $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.uncertainty

paste -d'\t' $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv" tmpfile.breakpoint.chr.end tmpfile.breakpoint.TE tmpfile.uncertainty > $SNa"_S"$SNo"_L"$LNo$"_out11_combined_results.tsv" &&

rm $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads.tsv"
rm tmpfile.*

echo " --> all done at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo".log"
