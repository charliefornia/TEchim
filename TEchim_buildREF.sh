#!/bin/bash

################################################################################
# TITLE: TEchim - build genomes
# VERSION: 0.1.1 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 13/05/2019 (dd/mm/yyyy)
# DESCRIPTION: tbd
################################################################################

################################################################################
# REQUIREMENTS:
# - bedtools
################################################################################

REF=~/Dropbox/CloudDesktop/REF_cloud
REFbase=dmel625_v04

cd $REF

################################################################################
# BUILD _GENES.bed
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$10); if ($3 == "gene") {print $1"\t"$4-1"\t"$5"\t"$10"\t"$6"\t"$7}}' | bedtools sort -i - > $REFbase"_GENES.bed"

################################################################################
# BUILD _SPLICE_DONORS.bed
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$10); if ($3 == "exon") {if ($7 == "+") {print $1"\t"sqrt(($5-11)^2)"\t"$5+10"\t"$10"|"$5"\t"$6"\t"$7} else if ($7 == "-") {print $1"\t"sqrt(($4-11)^2)"\t"$4+10"\t"$10"|"$4"\t"$6"\t"$7}}}' | bedtools sort -i - > $REFbase"_SPLICE_DONORS.bed"

################################################################################
# BUILD _SPLICE_ACCEPTORS.bed
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$10); if ($3 == "exon") {if ($7 == "-") {print $1"\t"sqrt(($5-11)^2)"\t"$5+10"\t"$10"|"$5"\t"$6"\t"$7} else if ($7 == "+") {print $1"\t"sqrt(($4-11)^2)"\t"$4+10"\t"$10"|"$4"\t"$6"\t"$7}}}' | bedtools sort -i - > $REFbase"_SPLICE_ACCEPTORS.bed"

################################################################################
# BUILD _INTRONS.gtf file - here, an intron is defined as any region between two
# exons of the same transcript.
# extract all lines that describe an "exon" in input gtf
cat $REFbase".gtf" | awk 'BEGIN {OFS = "\t"} {if ($3 == "exon") {print$0}}' > "tmp."$REFbase".exons.gtf"

# extract a (unique) list of all transcript names in input gtf
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$14); if ($3 == "exon") {print$14}}' | awk '!seen[$0]++' > "tmp."$REFbase".exons.transcripts.gtf"

# make sure that file does not yet exist in the current directory
rm -f $REFbase".output.introns.gtf"

# loop through each transcript in gtf file
while read line
do
	# if clause distinguishes between transcripts on the positive- and negative strand
	if [[ $(fgrep $line "tmp."$REFbase".exons.gtf" | head -n1 | awk '{print$7}') == "+" ]]
	then
		# on the positive strand, the start of the first intron is the end of the first exon minus one.
		# note that the last line of exons is deleted. this is because there's always one intron less than there are exons.
		LC_ALL=C fgrep $line "tmp."$REFbase".exons.gtf" | awk '{print $1"\t"$2"\tintron\t"$5+1}' | sed '$d' > "tmp."$REFbase".cola.tsv"
		# on the positive strand, the end of introns are the start of exons minus one.
		# note that here the first line is removed, not the last one. as a consequence, the values are shifted up one line.
		LC_ALL=C fgrep $line "tmp."$REFbase".exons.gtf" | awk 'NR>1 {print$4-1}' > "tmp."$REFbase".colb.tsv"
		# the remaining information about gene- and transcript name and id are added.
		LC_ALL=C fgrep $line "tmp."$REFbase".exons.gtf" | awk '{print $6"\t"$7"\t"$8"\t"$9"; "$10"; "$11"; "$12"; "$13"; "$14"; "$15"; "$16}' | sed '$d'  > "tmp."$REFbase".colc.tsv"
		paste "tmp."$REFbase".cola.tsv" "tmp."$REFbase".colb.tsv" "tmp."$REFbase".colc.tsv" >> $REFbase"_INTRONS.gtf"
		rm "tmp."$REFbase".col"*
	elif [[ $(fgrep $line "tmp."$REFbase".exons.gtf" | head -n1 | awk '{print$7}') == "-" ]]
	then
		# on the negative strand, the entries are stacked by removing the first line of the beginning and end of the entries.
		LC_ALL=C fgrep $line "tmp."$REFbase".exons.gtf" | awk 'NR>1 {print $1"\t"$2"\tintron\t"$5+1}' > "tmp."$REFbase".cola.tsv"
		# on the negative strand, the beginning of introns are the start of the exon minus one.
		LC_ALL=C fgrep $line "tmp."$REFbase".exons.gtf" | awk '{print$4-1}' | sed '$d' > "tmp."$REFbase".colb.tsv"
		# on the negative strand, the end of exons are the beginning of the trailing exon plus one.
		LC_ALL=C fgrep $line "tmp."$REFbase".exons.gtf" | awk 'NR>1 {print $6"\t"$7"\t"$8"\t"$9"; "$10"; "$11"; "$12"; "$13"; "$14"; "$15"; "$16}'   > "tmp."$REFbase".colc.tsv"
		paste "tmp."$REFbase".cola.tsv" "tmp."$REFbase".colb.tsv" "tmp."$REFbase".colc.tsv" >> $REFbase"_INTRONS.gtf"
		rm "tmp."$REFbase".col"*
	else
		echo "N/A\t"$line >> $REFbase"_INTRONS.gtf"
	fi
done < "tmp."$REFbase".exons.transcripts.gtf" &&\
rm "tmp."$REFbase*

################################################################################
# BUILD _FEATURES.bed
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "5UTR") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' > "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "CDS") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "stop_codon") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "3UTR") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "tRNA") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "pseudogene") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "snoRNA") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "ncRNA") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "snRNA") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "rRNA") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); if ($3 == "miRNA") {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> "tmp."$REFbase"_FEATURES.bed"
bedtools sort -i "tmp."$REFbase"_FEATURES.bed" > $REFbase"_FEATURES.bed"
bedtools sort -i $REFbase"_INTRONS.gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$0); {print $1"\t"$4-1"\t"$5"\t"$3"|"$10"|"$14"\t"$6"\t"$7}}' >> $REFbase"_FEATURES.bed"



