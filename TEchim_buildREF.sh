#!/bin/bash

################################################################################
# TITLE: TEchim - build genomes
# VERSION: 0.1.2 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 22/05/2019 (dd/mm/yyyy)
# DESCRIPTION: Run this once to create necesseary support files
################################################################################

################################################################################
# REQUIREMENTS:
# - bedtools
# - RepeatMasker
# - STAR (https://github.com/alexdobin/STAR)
# - blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
################################################################################

################################################################################
################################################################################
REFpath=/PATH/TO/REF/
REFbase=dmel625
TElist=TEs.fa
nc=10
# REFpath directory must contain:
#     - REFbase.fa (reference genome)
#     - REFbase.gtf
#     - TEs.fa (TE consensus sequences in fasta format)
################################################################################
################################################################################

# change to REF directory
cd $REFpath

# create TE.fa.bed file
awk '{if($1 !~ ">") {a=length($1); print a}}' $TElist > tmp.TElengths
awk '{if($1 ~ ">") {gsub(/>/,""); print $1}}' $TElist > tmp.TEnames
paste tmp.TElengths tmp.TEnames > tmp.TEcomb
awk '{print "TEchr_"$2"\t0\t"$1"\t"$2"\t.\t+"}' tmp.TEcomb > $TElist".bed"
awk '{print "TEchr_"$2"\t0\t"$1"\tNEG_"$2"\t.\t-"}' tmp.TEcomb >> $TElist".bed"
rm tmp.TE*

# use repeatmasker to mask any sequence in the reference genome that looks like a transposon
RepeatMasker -lib $TElist -no_is -nolow -s -pa $nc $REFbase".fa" 
rm $REFbase".fa.cat.gz"
rm $REFbase".fa.ori.out"
rm $REFbase".fa.out"
rm $REFbase".fa.tbl"
mv $REFbase".fa.masked" $REFbase".clean.noTEs.fa"

# create FASTA where each TE has chromosome name ">TEchr_..."
awk '{if ($1 ~ ">") {gsub(/>/,""); print ">TEchr_"$1"\t"$2"\t"$3} else {print $0}}' $TElist > $REFbase".clean.onlyTEs.fa"

# combine TE-cleaned reference genome with TE-fasta file
cat $REFbase".clean.noTEs.fa" $REFbase".clean.onlyTEs.fa" > $REFbase".clean.fa"

# generate STAR genome index
mkdir "STAR_"$REFbase
STAR --runMode genomeGenerate --genomeFastaFiles $REFbase".clean.fa" --genomeDir "./STAR_"$REFbase --runThreadN $nc
mv Log.out "STAR_"$REFbase"/."

# generate BLAST databases
makeblastdb -dbtype nucl -in $REFbase".clean.onlyTEs.fa"
makeblastdb -dbtype nucl -in $REFbase".clean.noTEs.fa"

# BUILD _GENES.bed
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$10); if ($3 == "gene") {print $1"\t"$4-1"\t"$5"\t"$10"\t"$6"\t"$7}}' | bedtools sort -i - > $REFbase"_GENES.bed"

# BUILD _EXONS.gtf
cat $REFbase".gtf" | awk '{if ($3 == "exon") {print $0}}' > $REFbase"_EXONS.gtf"

# BUILD _SPLICE_DONORS.bed
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$10); if ($3 == "exon") {if ($7 == "+") {print $1"\t"sqrt(($5-11)^2)"\t"$5+10"\t"$10"|"$5"\t"$6"\t"$7} else if ($7 == "-") {print $1"\t"sqrt(($4-11)^2)"\t"$4+10"\t"$10"|"$4"\t"$6"\t"$7}}}' | bedtools sort -i - > $REFbase"_SPLICE_DONORS.bed"

# BUILD _SPLICE_ACCEPTORS.bed
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$10); if ($3 == "exon") {if ($7 == "-") {print $1"\t"sqrt(($5-11)^2)"\t"$5+10"\t"$10"|"$5"\t"$6"\t"$7} else if ($7 == "+") {print $1"\t"sqrt(($4-11)^2)"\t"$4+10"\t"$10"|"$4"\t"$6"\t"$7}}}' | bedtools sort -i - > $REFbase"_SPLICE_ACCEPTORS.bed"

# BUILD _INTRONS.gtf file - here, an intron is defined as any region between two
# exons of the same transcript.
# extract all lines that describe an "exon" in input gtf
cat $REFbase".gtf" | awk 'BEGIN {OFS = "\t"} {if ($3 == "exon") {print$0}}' > "tmp."$REFbase".exons.gtf"
# extract a (unique) list of all transcript names in input gtf
cat $REFbase".gtf" | tr ';' '\t' | awk 'BEGIN {OFS = "\t"} {gsub(/\"/,"",$14); if ($3 == "exon") {print$14}}' | awk '!seen[$0]++' > "tmp."$REFbase".exons.transcripts.gtf"
# make sure that file does not yet exist in the current directory
rm -f $REFbase"_INTRONS.gtf"
# loop through each transcript in gtf file
while read line
do
	LC_ALL=C fgrep $line "tmp."$REFbase".exons.gtf" > "tmp."$REFbase".input.list"
	# if clause distinguishes between transcripts on the positive- and negative strand
	if [[ $(cat "tmp."$REFbase".input.list" | head -n1 | awk '{print$7}') == "+" ]]
	then
		# on the positive strand, the start of the first intron is the end of the first exon minus one.
		# note that the last line of exons is deleted. this is because there's always one intron less than there are exons.
		cat "tmp."$REFbase".input.list" | awk '{print $1"\t"$2"\tintron\t"$5+1}' | sed '$d' > "tmp."$REFbase".cola.tsv"
		# on the positive strand, the end of introns are the start of exons minus one.
		# note that here the first line is removed, not the last one. as a consequence, the values are shifted up one line.
		cat "tmp."$REFbase".input.list" | awk 'NR>1 {print$4-1}' > "tmp."$REFbase".colb.tsv"
		# the remaining information about gene- and transcript name and id are added.
		cat "tmp."$REFbase".input.list" | awk '{print $6"\t"$7"\t"$8"\t"$9"; "$10" "$11"; "$12" "$13"; "$14" "$15"; "$16}' | sed '$d'  > "tmp."$REFbase".colc.tsv"
		paste "tmp."$REFbase".cola.tsv" "tmp."$REFbase".colb.tsv" "tmp."$REFbase".colc.tsv" >> $REFbase"_INTRONS.gtf"
		rm "tmp."$REFbase".col"*
	elif [[ $(fgrep $line "tmp."$REFbase".exons.gtf" | head -n1 | awk '{print$7}') == "-" ]]
	then
		# on the negative strand, the entries are stacked by removing the first line of the beginning and end of the entries.
		cat "tmp."$REFbase".input.list" | awk 'NR>1 {print $1"\t"$2"\tintron\t"$5+1}' > "tmp."$REFbase".cola.tsv"
		# on the negative strand, the beginning of introns are the start of the exon minus one.
		cat "tmp."$REFbase".input.list" | awk '{print$4-1}' | sed '$d' > "tmp."$REFbase".colb.tsv"
		# on the negative strand, the end of exons are the beginning of the trailing exon plus one.
		cat "tmp."$REFbase".input.list" | awk 'NR>1 {print $6"\t"$7"\t"$8"\t"$9"; "$10" "$11"; "$12" "$13"; "$14" "$15"; "$16}'   > "tmp."$REFbase".colc.tsv"
		paste "tmp."$REFbase".cola.tsv" "tmp."$REFbase".colb.tsv" "tmp."$REFbase".colc.tsv" >> $REFbase"_INTRONS.gtf"
		rm "tmp."$REFbase".col"*
	else
		echo "N/A\t"$line >> $REFbase"_INTRONS.gtf"
	fi
done < "tmp."$REFbase".exons.transcripts.gtf"

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
rm "tmp."$REFbase*


# GENERATE support files for IGE
tmp_TEmin=$(grep -v ">" $TElist | awk '{l=length($1); print l}' | sort -n | head -n1)
cat $REFbase".gtf" | awk -v TEmin="$tmp_TEmin" '{if ($3 == "CDS" && $5-$4 > TEmin) {print $1"\t"$4-1"\t"$5"\t"$10"@"$14"\t.\t"$7}}' > "tmp."$REFbase".filtered_CDS.tsv"
bedtools getfasta -fi $REFbase.fa -bed "tmp."$REFbase".filtered_CDS.tsv" -name > "tmp."$REFbase".filtered_CDS.fasta"
makeblastdb -dbtype nucl -in $REFbase".fa"
blastn -query "tmp."$REFbase".filtered_CDS.fasta" -outfmt "10 qseqid" -db $REFbase".fa" | uniq -c | awk '{if($1 = "1") print $2}' > "tmp."$REFbase".use_these_CDS.tsv"
awk '{OFS="\t"} {gsub(/@|:/,"\t"); gsub(/"|;/,""); gsub(/-/,"\t",$5); print $4"\t"$5-1"\t"$6"\t"$1"@"$2"@"$4"-"$5":"$6"\t.\t"$3}' "tmp."$REFbase".use_these_CDS.tsv" > $REFbase".CDS_for_IGE.bed"
rm "tmp."$REFbase*
