#!/bin/bash

################################################################################
# TITLE: TEchim - PART 1
# VERSION: 0.3.0 (dev) - BRANCHED OUT AT THE MOMENT
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 07/02/2020 (dd/mm/yyyy)
# DESCRIPTION: This tool converts paired-end  OR single reads to in-silico
# pairs, which are then screened for pairs where one mate maps in the genome and
# the other on a transposon. The contigs of these read pairs are then used to 
# determine the precise breakpoint in the genome and the TE.
################################################################################

################################################################################
# REQUIREMENTS:
# - flash (https://ccb.jhu.edu/software/FLASH/)
# - STAR (https://github.com/alexdobin/STAR)
# - blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
# - samtools
################################################################################

################################################################################
################################################################################
# set parameters
wd=$(pwd)						# working directory
SNa=NAME_OF_EXP					# sample name
FASTQ1=READS_1.fastq.gz			# input FASTQ file 1 (FULL PATH)
FASTQ2=READS_2.fastq.gz			# input FASTQ file 2 (FULL PATH)
SNo=1							# sample number (use integer from 1-n)
LNo=1							# sequencing lane number (use integer from 1-n)
stranded=2						# strandedness: (1) FASTQ1 is mRNA strand
								#				(2) FASTQ2 is mRNA strand (default)
								#				(0) reads are unstranded
REFpath=/PATH/TO/REF/			# same path as in _buildREF
nc=1							# number of cores (default: 1)
fastalength=60					# length of in-silico FASTA - must be < 1/2 of
								# input read-length (default: 60)
################################################################################
################################################################################
# functions:

write_vars()
{
	# assign REFbase parameter
	if [ -e $REFpath"REFERENCE_basename" ]; then
		REFbase=$(cat $REFpath"REFERENCE_basename")
	else
		echo " #### ERROR: reference path is corrupt - no file named REFERENCE_basename" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
		exit
	fi
	if [ ! -e $wd"/."$SNa"_strandedness" ]; then echo $stranded > $wd"/."$SNa"_strandedness"; fi
	if [ ! -e $wd"/."$SNa"_samplename" ]; then echo $SNa > $wd"/."$SNa"_samplename"; fi
	if [ ! -e $wd"/."$SNa"_refpath" ]; then echo $REFpath > $wd"/."$SNa"_refpath"; fi
	if [ ! -e $wd"/."$SNa"_fastalength" ]; then echo $fastalength > $wd"/."$SNa"_fastalength"; fi
}

write_logfile()
{
	logname=$(date | awk '{gsub(/\:/,"-",$5); print $4$3$2"_"$5}')
	echo "====================" > $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "|| TEchim - PART1 || " >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "====================" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "Parameters:" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "Working directory:" "$wd/$SNa"_S"$SNo"_L"$LNo" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "FASTQ _1:" "$FASTQ1" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "FASTQ _2:" "$FASTQ2" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " => $(zcat -f < $FASTQ1 | wc -l | awk '{print $1/4}') reads" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "Sample name:" "$SNa" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "Sample number:" "$SNo" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "Sequencing lane number:" "$SNo" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "Reference files:" "$REFpath$REFbase" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "Length of in-silico reads:" "$fastalength""nt">> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "--------------------------------" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " --> starting at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "--------------------------------" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
}

merge_reads()
{
	echo " --> start merging at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	case $stranded in
		1|0)
			FQ1=$1
			FQ2=$2
			echo " ------ the first strand is the mRNA strand (unless input was unstranded)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
			;;
		2)
			FQ1=$2
			FQ2=$1
			echo " ------ the second strand is the mRNA strand" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
			;;
		*)
			echo " #### ERROR: incorrect strandedness: $stranded at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
			exit
			;;
	esac
	flash $FQ1 $FQ2 \
		-z \
		-x 0.15 \
		-M 170 \
		-o $SNa"_S"$SNo"_L"$LNo$"_out1" \
		-t $nc \
		-q
	# concatenate reads that were successfully merged with unmerged file _1.fastq
	# the unmerged file _1.fastq will contain reads that are longer than 2 x $fastalength,
	# which can also be used to generate in-silico paired-end reads
	cat $SNa"_S"$SNo"_L"$LNo$"_out1.extendedFrags.fastq.gz" $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_1.fastq.gz" > $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz"
	MaxFragLength=$(awk '{print $1}' $SNa"_S"$SNo"_L"$LNo$"_out1.hist" | tail -n1)
	echo $MaxFragLength >> $wd"/."$SNa"_maxfraglength"
	rm $SNa"_S"$SNo"_L"$LNo$"_out1.extendedFrags.fastq.gz"
	rm $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_1.fastq.gz"
	rm $SNa"_S"$SNo"_L"$LNo$"_out1.hist"
	rm $SNa"_S"$SNo"_L"$LNo$"_out1.histogram"
	echo " <-- done with merging at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "--------------------------------" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
}

create_fasta()
{ 
	echo " --> start generating cropped FASTA at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	# for the following awk commands, the reads will be filtered. only reads
	# containing at least 2 x $fastalength are kept. this is to avoid overlap when taking $fastalength
	# nucleotides from each end. 
	# first, split up fasta into 4 separate files, one for each line. (a) header;
	# (b) sequence; (c) +; (d) quality score.
	# reads used here: merged (_m) _1.fastq.gz (_1) -> _m1
	awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print a}}' <(gzip -dc $1) > a_m1.1
	awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print b}}' <(gzip -dc $1) > b_m1.2
	awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print c}}' <(gzip -dc $1) > c_m1.2
	awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print d}}' <(gzip -dc $1) > d_m1.2
	# add :A to readname of _m1 (This is to not confuse with otherwise
	# duplicate read names of _2)
	sed 's/ /:A /g' a_m1.1 > a_m1.2
	if [ -n $2 ]; then
		echo " ------ two FASTA files were used" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
		# reads used here: _2.fastq.gz reads (2)
		awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print a}}' <(gzip -dc $2) > a_2.1
		awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print b}}' <(gzip -dc $2) > b_2.1
		awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print c}}' <(gzip -dc $2) > c_2.2
		awk -v flength="$fastalength" 'BEGIN {OFS = "\n"} {a = $0 ; getline b ; getline c ; getline d ; if (length(b) >= (flength*2)) {print d}}' <(gzip -dc $2) > d_2.1
		# add :B to readname of _2 (This is to not confuse with otherwise
		# duplicate read names of _m1)
		sed 's/ /:B /g' a_2.1 > a_2.2
		rev b_2.1 | grep '^[ATCGN]' /dev/stdin | tr ATCGN TAGCN > b_2.2
		rev d_2.1 > d_2.2
	
		cat a_m1.2 a_2.2 > a_m12.1
		cat b_m1.2 b_2.2 > b_m12.1
		cat c_m1.2 c_2.2 > c_m12.1
		cat b_m1.2 d_2.2 > d_m12.1
	else
		cat a_m1.2 > a_m12.1
		cat b_m1.2 > b_m12.1
		cat c_m1.2 > c_m12.1
		cat b_m1.2 > d_m12.1	
	fi
	# crop sequences and quality scores to $fastalength - for out3_1.fasta
	sed -E "s/(.{$fastalength}).*/\1/" b_m12.1 > b_m12.2
	sed -E "s/(.{$fastalength}).*/\1/" d_m12.1 > d_m12.2
	# crop sequences ($fastalength) and generate reverse-complement - for out3_2.fasta
	rev b_m12.1 | sed -E "s/(.{$fastalength}).*/\1/" | grep '^[ATCGN]' /dev/stdin | tr ATCGN TAGCN > b_m12.3
	# crop quality scores ($fastalength) and invert - for out3_2.fasta
	rev d_m12.1 | sed -E "s/(.{$fastalength}).*/\1/" > d_m12.3
	# Create a in-silico FASTA and lookup table for long reads. in-silico_1.fasta and the long read  of LOOKUP
	# is the ACTUAL nucleotide sequence of the mRNA molecule 
	sed 's/ .*//' a_m12.1 > a_m12.2
	# combine in-silico components
	paste -d'\n' a_m12.1 b_m12.2 c_m12.1 d_m12.2 > $SNa"_S"$SNo"_L"$LNo$"_out3_1.fasta"
	paste -d'\n' a_m12.1 b_m12.3 c_m12.1 d_m12.3 > $SNa"_S"$SNo"_L"$LNo$"_out3_2.fasta"
	# create look-up table and remove "@" at the beginning of the readname
	paste -d'\t' a_m12.2 b_m12.1 | sed 's/^@\(.*\)/\1/' > $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.tsv"
	# create FASTA file from LOOKUP table
	awk '{print ">"$1"\n"$2}' $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.tsv" > $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.fa"
	# LANG=en_EN is a bug-fix to make sort compatible with join
	LANG=en_EN sort -k 1,1 $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.tsv" | gzip > $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.sorted.tsv.gz" && rm $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.tsv"
	rm a_*
	rm b_*
	rm c_*
	rm d_*
	echo " ------ in-silico FASTA contain $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out3_1.fasta" | awk '{print $1/4}') reads" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " <-- done with generating FASTA at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "--------------------------------" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
}

align_and_filter()
{
	echo " --> start mapping at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	#### STREAM 1 ####
	# run STAR in chimera-mode
	STAR --runThreadN $nc --genomeDir $REFpath"STAR_"$REFbase --readFilesIn $1 $2 --chimSegmentMin 20 --chimOutType WithinBAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $SNa"_S"$SNo"_L"$LNo$"_out4_" --clip3pAdapterSeq AAAAAAAA
	mkdir $SNa"_S"$SNo"_L"$LNo"_STAR"
	mv $SNa"_S"$SNo"_L"$LNo$"_out4"* $SNa"_S"$SNo"_L"$LNo"_STAR"/.
	# extract only hits that cross TE-GENE breakpoints. the awk commands remove
	# (1) TE-TE reads and (2)&(3) reads that span the TE and it's LTR
	samtools view $SNa"_S"$SNo"_L"$LNo"_STAR"/$SNa"_S"$SNo"_L"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | grep TEchr_ | awk '($7 != "=") && ($3 !~ "TEchr_" || $7 !~ "TEchr_")' > $SNa"_S"$SNo"_L"$LNo$"_out5_STREAM1_TExGenes.sam"
	if [ ! -s $SNa"_S"$SNo"_L"$LNo$"_out5_STREAM1_TExGenes.sam" ]; then
		echo " #### ERROR: file does not have any TE-GENE brakpoint spanning reads!" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
		echo " --> exited at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
		exit
	fi
	#### STREAM 2 ####
	# The output of this stream is based on the STAR alignment of in-silico paired-end reads.
	# The hits will only be used for reads where the genome-section is successfully mapped with BLAST (further downstream), but the transposon section is NOT.
	# Using the maximum fragment length, and information about whether the gene- and te- reads are mapped to the (+)ive or (-)ive strand (all contained in SAM-flag) are used for output information.
	# Select pairs where both reads map to positive strand
	samtools view -f 65 -F 48 $SNa"_S"$SNo"_L"$LNo"_STAR"/$SNa"_S"$SNo"_L"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4,$4+mfl,$1"|"$7"|minus|GENE-TE|S"s"|L"l"|"$8,".","+",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8,$8+mfl,$1"|"$3"|plus|TE-GENE|S"s"|L"l"|"$4,".","-",$1}}' | sed 's/TEchr_//g' > $SNa"_S"$SNo"_L"$LNo$"_out5b_STREAM2_FORout10c.bed"
	# Select pairs where read1 maps to positive, and read2 to negative
	samtools view -f 97 -F 16 $SNa"_S"$SNo"_L"$LNo"_STAR"/$SNa"_S"$SNo"_L"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4,$4+mfl,$1"|"$7"|plus|GENE-TE|S"s"|L"l"|"$8,".","+",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8-mfl,$8,$1"|"$3"|plus|TE-GENE|S"s"|L"l"|"$4,".","+",$1}}' | sed 's/TEchr_//g' >> $SNa"_S"$SNo"_L"$LNo$"_out5b_STREAM2_FORout10c.bed"
	# Select pairs where read1 maps to negative, and read2 to positive
	samtools view -f 81 -F 32 $SNa"_S"$SNo"_L"$LNo"_STAR"/$SNa"_S"$SNo"_L"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4-mfl,$4,$1"|"$7"|minus|GENE-TE|S"s"|L"l"|"$8,".","-",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8,$8+mfl,$1"|"$3"|minus|TE-GENE|S"s"|L"l"|"$4,".","-",$1}}' | sed 's/TEchr_//g' >> $SNa"_S"$SNo"_L"$LNo$"_out5b_STREAM2_FORout10c.bed"
	# Select pairs where both reads map to negative strand
	samtools view -f 113 $SNa"_S"$SNo"_L"$LNo"_STAR"/$SNa"_S"$SNo"_L"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4-mfl,$4,$1"|"$7"|plus|GENE-TE|S"s"|L"l"|"$8,".","-",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8-mfl,$8,$1"|"$3"|minus|TE-GENE|S"s"|L"l"|"$4,".","+",$1}}' | sed 's/TEchr_//g' >> $SNa"_S"$SNo"_L"$LNo$"_out5b_STREAM2_FORout10c.bed"
	bedtools sort -i $SNa"_S"$SNo"_L"$LNo$"_out5b_STREAM2_FORout10c.bed" > $SNa"_S"$SNo"_L"$LNo$"_out10b_STREAM2_additional_chimera.bed"
	rm $SNa"_S"$SNo"_L"$LNo$"_out5b_STREAM2_FORout10c.bed"
	rm $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.fa"
	echo " ------ sample contains $(awk '{print $1}' $SNa"_S"$SNo"_L"$LNo$"_out5_STREAM1_TExGenes.sam" | sort | uniq | wc -l | awk '{print $1}') unique reads that span gene-TE breakpoint." >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " ------ sample contains $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out10b_STREAM2_additional_chimera.bed" | awk '{print $1}') unique reads that were picked up for STREAM2 - these will still be filtered, STREAM1 has priority." >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " ------ sample contains $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out10c_STREAM3_additional_breakpoints.unfiltered.bed" | awk '{print $1}') unique reads that were picked up for STREAM3 - these will still be filtered, STREAM1 and 2 have priority." >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " <-- done with mapping at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "--------------------------------" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
}

blast_on_longreads ()
{
	echo " --> start BLAST at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	# extract readnames and remove duplicates
	awk '{print $1}' $1 | sort | uniq > $SNa"_S"$SNo"_L"$LNo$"_out6_TExGENES_readnames.txt"
	# combine readnames with long sequences stored in the lookup file
	# LANG=en_EN is a bug-fix to make sort compatible with join
	LANG=en_EN join -1 1 -2 1 <(sort $SNa"_S"$SNo"_L"$LNo$"_out6_TExGENES_readnames.txt") <(zcat < $2 ) > $SNa"_S"$SNo"_L"$LNo$"_out7_TExGENES_longreads.tsv" && rm $SNa"_S"$SNo"_L"$LNo$"_out6_TExGENES_readnames.txt"
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
		echo "${sequence}" | blastn -db $REFpath$REFbase".clean.noTEs.fa" -outfmt "6 qstart qend sseqid sstart send sstrand qlen" -num_alignments 1 -num_threads $nc | head -n 1 > var3		
		if ! [ -s var3 ]; then echo "noGENEfound" > var3; fi
		echo "${sequence}" | blastn -db $REFpath$REFbase".clean.onlyTEs.fa" -outfmt "6 qstart qend sseqid sstart send slen sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var4
		if ! [ -s var4 ]; then echo "noTEfound" > var4; fi
		paste var1 var2 var3 var4 >> $SNa"_S"$SNo"_L"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
		rm var*
	done < $SNa"_S"$SNo"_L"$LNo$"_out7_TExGENES_longreads.tsv"
	# remove reads that did not give BLAST result for genomic location
	grep -v noGENEfound $SNa"_S"$SNo"_L"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv" > $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv"
	# extract reads for which TE has been identified
	grep -v noTEfound $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv" | awk 'BEGIN {OFS = "\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16}' > $SNa"_S"$SNo"_L"$LNo$"_out10a_STREAM1_TExGENES_blastedreads.tsv"
	# extract reads where NO TE was found - these are filtered (at least $fastalength/2 should NOT map to genome) and STREAM2 information is added
	# remove reads where less than half the $fastalength is left unmapped when BLASTed to no_TE_genome	
	grep noTEfound $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv" | awk -v flength="$fastalength" '{if ($9-$4 > flength/2 || $3 > flength/2) {print $0} }' > $SNa"_S"$SNo"_L"$LNo$"_out9b_fromSTREAM1forSTREAM2_TExGENES_findTE.tsv"
	join -1 1 -2 7 <(sort $SNa"_S"$SNo"_L"$LNo$"_out9b_fromSTREAM1forSTREAM2_TExGENES_findTE.tsv") <(sort -k 7 $SNa"_S"$SNo"_L"$LNo$"_out10b_STREAM2_additional_chimera.bed") | awk 'BEGIN {OFS = "\t"} {if ($8=="plus") {if ($14 ~ "[\|]TE-GENE[\|]") {print $5,$6-1,$6,$14,".","+",$3,$4,$9} else if ($14 ~ "[\|]GENE-TE[\|]") {print $5,$7-1,$7,$14,".","+",$3,$4,$9}} else if ($8=="minus") {if ($14 ~ "[\|]TE-GENE[\|]") {print $5,$6-1,$6,$14,".","-",$3,$4,$9} else if ($14 ~ "[\|]GENE-TE[\|]") {print $5,$7-1,$7,$14,".","-",$3,$4,$9}}}' | bedtools sort -i - > $SNa"_S"$SNo"_L"$LNo$"_out11b_STREAM2_additional_chimera.filtered.bed"
	tr "|" "\t" < $SNa"_S"$SNo"_L"$LNo$"_out11b_STREAM2_additional_chimera.filtered.bed" | awk -v flength="$fastalength" 'BEGIN {OFS="\t"} {if ($6 == "plus") {if ($7 == "TE-GENE") {testart=$10+flength; teend=$10+$13} else if ($7 == "GENE-TE") {testart=$10+flength+$14-$15; teend=$10}} else if ($6 == "minus") {if ($7 == "TE-GENE") {testart=$10+flength-$13; teend=$10} else if ($7 == "GENE-TE") {testart=$10+flength; teend=$10+$15-$14}} {print $1,$2,$3,$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"testart"-"teend"|0",$11,$12}}' > $SNa"_S"$SNo"_L"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed"
	echo " ------ BLAST results:" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " ------ $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out7_TExGENES_longreads.tsv" | awk '{print $1}') input reads" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " ------ $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv" | awk '{print $1}') reads with a genome location" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " ------ |> of the reads with genome location: $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out10a_STREAM1_TExGENES_blastedreads.tsv" | awk '{print $1}') reads with a genome- and transposon location" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo " ------ |> of the reads with genome location: $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed" | awk '{print $1}') reads from STREAM2 will be added to STREAM1" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out7_TExGENES_longreads.tsv"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out9b_fromSTREAM1forSTREAM2_TExGENES_findTE.tsv"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out10b_STREAM2_additional_chimera.bed"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out11b_STREAM2_additional_chimera.filtered.bed"
	echo " <-- done with BLAST at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	echo "--------------------------------" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
}

create_summary_table ()
{
	echo " --> start creating summary table at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	# determine whether the section that maps to the genome is the 5' or the 3' end
	# of the mRNA section:
	# 5'-######|GENE|TE|######-3' => PART5
	# 5'-######|TE|GENE|######-3' => PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; if (a < b) {print "PART5"} else {print "PART3"}}' < $1 > tmpfile.genepart
	# determine the precise breakpoint on the chromosome. this depends on  whether
	# the chromosomal part is PART5 or PART3. BLAST places the read that matches $3 to $6, regardless of the orientation.
	awk 'BEGIN {OFS = "\t"} {a = $5; {print a}}' < $1 > tmpfile.chr
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d-1} else {print c-1}}' < $1 > tmpfile.breakpoint.chr.start
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d} else {print c}}' < $1 > tmpfile.breakpoint.chr.end
	# determine the precise breakpoint on the TE. this depends on  whether the TE
	# part is PART5 or PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $13 ; d = $12 ; if (a < b) {print d} else {print c}}' < $1 > tmpfile.breakpoint.TE
	# determine the overlap between the two mapped sections of the long read
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $4 ; d = $10 ; if (a < b) {print b-c-1} else {print a-d-1}}' < $1 > tmpfile.uncertainty
	if [[ $stranded = "0" ]]; then
		awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
			a = $11
			gsub(/TEchr_/,"",a)
			if ($8 == "plus") {if ($15=="plus") {b = "forward"} else {b = "reverse"}} else {if ($15 == "plus") { b = "reverse" } else { b = "forward" }}
			if ($3 < $9) {print $1"|"a"|"b"|GENE-TE|S"s"|L"l} else {print $1"|"a"|"b"|TE-GENE|S"s"|L"l}
			}' < $1 > tmpfile.readname
			awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "+"}}' < $1 > tmpfile.breakpoint.chr.strand
	else
		awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
		a = $11
		gsub(/TEchr_/,"",a)
		b = $15
		c = $3
		d = $9
		e = $1
		if (c < d) {print e"|"a"|"b"|GENE-TE|S"s"|L"l} else {print e"|"a"|"b"|TE-GENE|S"s"|L"l}
		}' < $1 > tmpfile.readname
		awk 'BEGIN {OFS = "\t"} {a = $8 ; if (a == "plus") {print "+"} else {print "-"}}' < $1 > tmpfile.breakpoint.chr.strand
	fi
	awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "."}}' < $1 > tmpfile.score
	paste -d'|' tmpfile.readname tmpfile.breakpoint.TE tmpfile.uncertainty > tmpfile.readname.extended
	paste -d'\t' tmpfile.chr tmpfile.breakpoint.chr.start tmpfile.breakpoint.chr.end tmpfile.readname.extended tmpfile.score tmpfile.breakpoint.chr.strand > $SNa"_S"$SNo"_L"$LNo$"_out11a_STREAM1_breakpoints.bed"
	# add results from STREAM 2
	cat $SNa"_S"$SNo"_L"$LNo$"_out11a_STREAM1_breakpoints.bed" $SNa"_S"$SNo"_L"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed" | bedtools sort -i - > $SNa"_S"$SNo"_L"$LNo$"_out13_breakpoints.bed"
	echo " ------ In total, $(wc -l $SNa"_S"$SNo"_L"$LNo$"_out13_breakpoints.bed" | awk '{print $1}') chimeric reads were found." >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	rm tmpfile.*
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out11a_STREAM1_breakpoints.bed"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed"
	rm -f $SNa"_S"$SNo"_L"$LNo$"_out12_STREAM1and2_breakpoints.bed"
}

################################################################################
################################################################################
# execute:

# change to wd
cd $wd
# change to folder for sample
if [ ! -d $SNa"_S"$SNo"_L"$LNo ]; then
	mkdir $SNa"_S"$SNo"_L"$LNo
fi
cd $SNa"_S"$SNo"_L"$LNo

write_vars
write_logfile
if [ -f "$FASTQ2" ]; then
	merge_reads $FASTQ1 $FASTQ2
	create_fasta $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz" $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_2.fastq.gz" \
		&& rm $SNa"_S"$SNo"_L"$LNo$"_out1.combined.fastq.gz" && rm $SNa"_S"$SNo"_L"$LNo$"_out1.notCombined_2.fastq.gz"
elif [ -f "$FASTQ1" ]; then
	create_fasta $FASTQ1
	MaxFragLength=$(awk '{ if (length($0) > max) max = length($0) } END { print max }' $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.fa")
	echo $MaxFragLength >> $wd"/."$SNa"_maxfraglength"
else
	echo " #### ERROR: At least one FASTQ input is required!" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
	exit
fi

align_and_filter $SNa"_S"$SNo"_L"$LNo$"_out3_1.fasta" $SNa"_S"$SNo"_L"$LNo$"_out3_2.fasta"
blast_on_longreads $SNa"_S"$SNo"_L"$LNo$"_out5_STREAM1_TExGenes.sam" $SNa"_S"$SNo"_L"$LNo$"_LOOKUP.sorted.tsv.gz" &&\
	rm $SNa"_S"$SNo"_L"$LNo$"_out5_STREAM1_TExGenes.sam"

create_summary_table $SNa"_S"$SNo"_L"$LNo$"_out10a_STREAM1_TExGENES_blastedreads.tsv" &&\
	rm $SNa"_S"$SNo"_L"$LNo$"_out10a_STREAM1_TExGENES_blastedreads.tsv"

echo " <-- all done at ... $(date)" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
echo "================================" >> $SNa"_S"$SNo"_L"$LNo"_PART1_"$logname".log"
