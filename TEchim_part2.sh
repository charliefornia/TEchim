#!/bin/bash

################################################################################
# TITLE: TEchim - part 2
# VERSION: 0.1.1 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 09/05/2019 (dd/mm/yyyy)
# DESCRIPTION: tbd
################################################################################

################################################################################
# REQUIREMENTS:
# - bedtools
################################################################################

# set parameters
wd=~
path_to_PART1_output=~/PATH/TO/PART1/OUTPUT/
SNa=NameOfSample
REF=~/PATH/TO/REF/
REFbase=dmel625_v04

# filters for "anchor" insertion
min_repl=2
min_reads=4

# change to wd
cd $wd

# create out10_combined.sorted.bed
cat $path_to_PART1_output$SNa*"/"$SNa*"_out10_breakpoints.bed" | bedtools sort -i - > $SNa"_in10_combined.sorted.bed"

# add a gene tag to each read (i.e. genomic location). this is done using a gene-only version of the gtf file.
bedtools intersect -wa -a $SNa"_in10_combined.sorted.bed" -b $REF$REFbase"_GENES.bed" -loj -s > $SNa"_out01_genetagged.tsv"

# separate | delimited field
sed -e $'s/|/\t/g' $SNa"_out01_genetagged.tsv" > $SNa"_out02_sepparated.tsv" &&\
	rm $SNa"_out01_genetagged.tsv"

# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
awk 'BEGIN {OFS = "\t"} {
	a = $5
	gsub(/_LTR/,"",a)
	c = substr($4, 1, length($4)-2)
	print $0"\t"a"\t"c
	}' < $SNa"_out02_sepparated.tsv" > $SNa"_out03pre_TEbase.tsv" &&\
		rm $SNa"_out02_sepparated.tsv"

# remove duplicate reads (where both :A and :B version of the same reads were picked up
awk '!seen[$21]++' $SNa"_out03pre_TEbase.tsv" > $SNa"_out03_TEbase.tsv" &&\
	rm $SNa"_out03pre_TEbase.tsv"

# create list of all TEs in dataset
cut -f20 $SNa"_out03_TEbase.tsv" | sort | uniq > $SNa"_out03a_uniqueTEs.tsv"

# to create collection file, first make sure this file does not yet exist
rm -f $SNa"_out06.maininsertions.tsv"
rm -f $SNa"_out12_OUTPUT.tsv"

# loop through all TEs in the dataset
while read TE
do
	# grep TEs from combined data set
	awk -v t="$TE" '{if ($20 == t) print $0;}' $SNa"_out03_TEbase.tsv" > "tmp."$SNa"_"$TE"_out04.tsv"
	
	###### SPLIT into TE-GENE and GENE-TE
	awk '{if ($7 == "TE-GENE") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv"
	awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv"
	
	# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
	# header: Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllBreakpoints(TE)
	bedtools merge -i "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_out06_TE-GENE.tsv"
	bedtools merge -i "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_out06_GENE-TE.tsv"
	
########### FOR ANCHOR
	# filter to keep only insertions that were detected in at least $min_repl separate biological replicates with at least $min_reads reads
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {print a}}' < "tmp."$SNa"_"$TE"_out06_TE-GENE.tsv" > "tmp."$SNa"_"$TE"_out07_TE-GENE.tsv"
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {print a}}' < "tmp."$SNa"_"$TE"_out06_GENE-TE.tsv" > "tmp."$SNa"_"$TE"_out07_GENE-TE.tsv"
	# combine TE-GENE and GENE-TE, then sort
	cat "tmp."$SNa"_"$TE"_out07_TE-GENE.tsv" "tmp."$SNa"_"$TE"_out07_GENE-TE.tsv" > "tmp."$SNa"_"$TE"_out08.tsv"
	bedtools sort -i "tmp."$SNa"_"$TE"_out08.tsv" > "tmp."$SNa"_"$TE"_out09.tsv"

########### FOR rest
	# filter to keep only insertions that were detected in at least $min_repl separate biological replicates with at least $min_reads reads
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {} else {print a}}' < "tmp."$SNa"_"$TE"_out06_TE-GENE.tsv" > "tmp."$SNa"_"$TE"_out07b_TE-GENE_REST.tsv"
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {} else {print a}}' < "tmp."$SNa"_"$TE"_out06_GENE-TE.tsv" > "tmp."$SNa"_"$TE"_out07b_GENE-TE_REST.tsv"
	# combine TE-GENE and GENE-TE, then sort
	cat "tmp."$SNa"_"$TE"_out07b_TE-GENE_REST.tsv" "tmp."$SNa"_"$TE"_out07b_GENE-TE_REST.tsv" > "tmp."$SNa"_"$TE"_out08b_REST.tsv"
	bedtools sort -i "tmp."$SNa"_"$TE"_out08b_REST.tsv" > "tmp."$SNa"_"$TE"_out09b_REST.tsv"


	# create list of all genes in ACHOR file
	# first, replace "," with tab, then sort, then pick unique strings, then remove ""
	cut -f12 "tmp."$SNa"_"$TE"_out09.tsv" | tr ',' '\n' | sort | uniq | sed '/\./d' > "tmp."$SNa"_"$TE"_out11_genes.tsv"
	
	# this list of genes is used to add non-anchor hits to ANCHOR hits, if they overlap with the same gene and TE
	while read GENE
		do
			# create output file that contains for every gene with at least one anchor insertion, all the insertions of that TE
			# ANCHOR hits have "yes" appended, others "NO"
			awk -v g="$GENE" '{if ($12 ~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$13"\t"$7"\t"$8"\t""YES";}' "tmp."$SNa"_"$TE"_out09.tsv" >> $SNa"_out12_OUTPUT.tsv"
			awk -v g="$GENE" '{if ($12 ~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$13"\t"$7"\t"$8"\t""NO";}' "tmp."$SNa"_"$TE"_out09b_REST.tsv" >> $SNa"_out12_OUTPUT.tsv"
		done < "tmp."$SNa"_"$TE"_out11_genes.tsv"
	
	rm "tmp."$SNa*
done < $SNa"_out03a_uniqueTEs.tsv"

rm $SNa"_out03a_uniqueTEs.tsv"
rm $SNa"_out03_TEbase.tsv"

# now, every line is one incident of a gene-TE chimera
while read line
do
	# convert line to bed format
	# header: Chr(genome)\tStart(genome)[this is the breakpoint minus 1]\tEnd(genome)[this is the breakpoint]\tTE|Gene|TE-GENEorGENE-TE\t"."\tStrand(genome)
	echo $line | awk '{ print $2"\t"$4-1"\t"$4"\t"$5"|"$1"|"$6"\t.\t"$3}' > "tmp."$SNa"_out13.bed"
	
	# $a is TE-GENE or GENE-TE	
	a=$(echo $line | awk '{print $6}')
	# $b is the gene
	b=$(echo $line | awk '{print $1}')
	
	# intersect with FEATURES file to get all ovelapping features (output is ;-separated)
	d=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b $REF$REFbase"_FEATURES.bed" -loj -wa -s | awk '{print $10}' | paste -sd ";" -)
	
	# check whether breakpoint on the genome overlaps with splice donor site (in the case of a GENE-TE fragment)
	# or with splice acceptor site (in the case of a TE-GENE fragment). the REF files _SPLICE_DONORS.bed and 
	# _SPLICE_ACCEPTORS.bed have a +/- 10nt, but the precise breakpoint is extracted here (which is the sequence
	# number of the last nucleotide that is still part of the exon for SPLICE_DONORS and the first exonic nt 
	# for SPLICE_ACCEPTORS.
	# for GENE-TE fragments, the exon is the splice donor
	if [ $a == "GENE-TE" ]
	then
		grep $b $REF$REFbase"_SPLICE_DONORS.bed" > "tmp."$SNa"_out14_ref_for_overlap.bed"
		c=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b "tmp."$SNa"_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
		rm "tmp."$SNa"_out14_ref_for_overlap.bed"
	# for TE-GENE fragments, the exon is the splice acceptor
	elif [ $a == "TE-GENE" ]
	then
		grep $b $REF$REFbase"_SPLICE_ACCEPTORS.bed" > "tmp."$SNa"_out14_ref_for_overlap.bed"
		c=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b "tmp."$SNa"_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
		rm "tmp."$SNa"_out14_ref_for_overlap.bed"
	else
		c="N/A"
	fi
	# the loop aboce created the unix variable $c. here, this variable is appended to the main file.
	echo $line | awk -v c="$c" -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"c"\t"d}'
	rm "tmp."$SNa"_out13.bed"
done < $SNa"_out12_OUTPUT.tsv" > $SNa"_out15.tsv" &&\
	rm $SNa"_out12_OUTPUT.tsv"

# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)
cat $SNa"_out15.tsv" | while read line; do echo $line | awk '{print $8}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| sed 's/,$//'; done > $SNa"_newcolb"
# for final output, get rid of original TE-breakpoint field
awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13;}' $SNa"_out15.tsv" > $SNa"_newcola" &&\
	rm $SNa"_out15.tsv"
# append pooled TE breakpoints to final output
paste $SNa"_newcola" $SNa"_newcolb" > $SNa"_chimericreads_final.tsv"
rm $SNa"_newcola" $SNa"_newcolb"
