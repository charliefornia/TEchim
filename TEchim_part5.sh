#!/bin/bash

################################################################################
# TITLE: TEchim - PART 5
# VERSION: 0.2.0 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 04/06/2019 (dd/mm/yyyy)
# DESCRIPTION: This script combines the IGE analysis from separate samples and
# lanes from PART4 and generates summary .tsv files with all necessary data
################################################################################

################################################################################
# REQUIREMENTS:
# - bedtools
# - samtools
################################################################################

################################################################################
################################################################################
# set parameters
wd=~
path_to_PART4_output=/PATH/TO/PART4/OUTPUT/
SNa=NAME_OF_EXP
REFpath=/PATH/TO/REF/
REFbase=dmel625
################################################################################
################################################################################
# functions:

process_PART1_output()
{
	echo " --> start processing PART1 output at ... $(date)" >> $SNa"_PART5_"$logname".log"
	# filter reads inside same gene
	awk '{if ($4 !~ $10) print $0}' $1 > $SNa"_IGEref_out01a_filtered_genetagged.tsv"
	# separate | delimited field
	sed -e $'s/|/\t/g' $SNa"_IGEref_out01a_filtered_genetagged.tsv" > $SNa"_IGEref_out02_sepparated.tsv" && rm $SNa"_IGEref_out01a_filtered_genetagged.tsv" && rm $SNa"_IGEref_out01_genetagged.tsv"
	# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
	# and create column with BASE readname (no :A or :B)
		awk 'BEGIN {OFS = "\t"} {
		a = $5
		gsub(/_LTR/,"",a)
		c = substr($4, 1, length($4)-2)
		print $0"\t"a"\t"c
		}' < $SNa"_IGEref_out02_sepparated.tsv" > $SNa"_IGEref_out03pre_TEbase.tsv" && rm $SNa"_IGEref_out02_sepparated.tsv"
	wc_predup=$(wc -l $SNa"_IGEref_out03pre_TEbase.tsv" | awk '{print $1}')
	# remove duplicate reads (where both :A and :B version of the same reads were picked up
	awk '!seen[$21]++' $SNa"_IGEref_out03pre_TEbase.tsv" > $SNa"_IGEref_out03_TEbase.tsv" && rm $SNa"_IGEref_out03pre_TEbase.tsv"
	echo " ---- $(wc -l $SNa"_IGEref_out03_TEbase.tsv" | awk '{print $1}') reads ($wc_predup before removing duplicates)" >> $SNa"_PART5_"$logname".log"
	echo " <-- done processing PART1 output at ... $(date)" >> $SNa"_PART5_"$logname".log"
}

combine_hits_of_each_TE()
{
	echo " --> start combining reads at ... $(date)" >> $SNa"_PART5_"$logname".log"
	# create list of all TEs in dataset
	cut -f20 $1 | sort | uniq > $SNa"_IGEref_out03a_uniqueTEs.tsv"
	# to create collection file, first make sure this file does not yet exist
	rm -f $SNa"_IGEref_out12_OUTPUT.tsv"
	# loop through all TEs in the dataset
	while read TE
	do
		# grep TEs from combined data set
		awk -v t="$TE" '{if ($20 == t) print $0;}' $1 > "tmp."$SNa"_"$TE"_IGEref_out04.tsv"
		if [ -s "tmp."$SNa"_"$TE"_IGEref_out04.tsv" ]; then	
			###### SPLIT into TE-GENE and GENE-TE
			awk '{if ($7 == "TE-GENE") print $0;}' "tmp."$SNa"_"$TE"_IGEref_out04.tsv" > "tmp."$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv"
			# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
			# header: Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllBreakpoints(TE)
			if [ -s "tmp."$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_IGEref_out06_TE-GENE.tsv"
			fi
			awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$SNa"_"$TE"_IGEref_out04.tsv" > "tmp."$SNa"_"$TE"_IGEref_out05_GENE-TE.tsv"
			if [ -s "tmp."$SNa"_"$TE"_IGEref_out06_GENE-TE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_IGEref_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_IGEref_out06_GENE-TE.tsv"
			fi
			###### COMBINE TE-GENE and GENE-TE
			cat "tmp."$SNa"_"$TE"_IGEref_out06_"*".tsv" > "tmp."$SNa"_"$TE"_IGEref_out08.tsv"
			bedtools sort -i "tmp."$SNa"_"$TE"_IGEref_out08.tsv" > "tmp."$SNa"_"$TE"_IGEref_out09.tsv"
			# create list of all genes
			# first, replace "," with tab, then sort, then pick unique strings, then remove ""
			cut -f12 "tmp."$SNa"_"$TE"_IGEref_out09.tsv" | tr ',' '\n' | sort | uniq | sed '/\./d' > "tmp."$SNa"_"$TE"_IGEref_out11_genes.tsv"
			# this list of genes is used to add non-anchor hits to ANCHOR hits, if they overlap with the same gene and TE
			while read GENE
			do
				# create output file that contains for every gene all the insertions of that TE
				awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$13"\t"$7"\t"$8"\t""mRNA";}' "tmp."$SNa"_"$TE"_IGEref_out09.tsv" >> $SNa"_IGEref_out12_OUTPUT.tsv"
			done < "tmp."$SNa"_"$TE"_IGEref_out11_genes.tsv"
			rm -f "tmp."$SNa"_"$TE*
		fi
	done < $SNa"_IGEref_out03a_uniqueTEs.tsv"
	# remove empty lines
	grep . $SNa"_IGEref_out12_OUTPUT.tsv" > $SNa"_IGEref_out12a_OUTPUT.tsv"
	rm $SNa"_IGEref_out03a_uniqueTEs.tsv"
	rm $SNa"_IGEref_out03_TEbase.tsv"
	echo " <-- done combining reads at ... $(date)" >> $SNa"_PART5_"$logname".log"
}

check_for_splicesites()
{
	echo " --> start checking for breaks at splice sites at ... $(date)" >> $SNa"_PART5_"$logname".log"
	# now, every line is one incident of a gene-TE chimera
	while read line
	do
		# convert line to bed format
		# header: Chr(genome)\tStart(genome)[this is the breakpoint minus 1]\tEnd(genome)[this is the breakpoint]\tTE|Gene|TE-GENEorGENE-TE\t"."\tStrand(genome)
		echo $line | awk '{ print $2"\t"$4-1"\t"$4"\t"$5"|"$1"|"$6"\t.\t"$3}' > "tmp."$SNa"_IGEref_out13.bed"
		# $a is TE-GENE or GENE-TE	
		a=$(echo $line | awk '{print $6}')
		# $b is the gene
		b=$(echo $line | awk '{print $1}')
		# intersect with FEATURES file to get all ovelapping features (output is ;-separated)
		d=$(bedtools intersect -a "tmp."$SNa"_IGEref_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa -s | awk '{print $10}' | paste -sd ";" -)
		# check whether breakpoint on the genome overlaps with splice donor site (in the case of a GENE-TE fragment)
		# or with splice acceptor site (in the case of a TE-GENE fragment). the REF files _SPLICE_DONORS.bed and 
		# _SPLICE_ACCEPTORS.bed have a +/- 10nt, but the precise breakpoint is extracted here (which is the sequence
		# number of the last nucleotide that is still part of the exon for SPLICE_DONORS and the first exonic nt 
		# for SPLICE_ACCEPTORS.
		# for GENE-TE fragments, the exon is the splice donor
		if [ $a == "GENE-TE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_DONORS.bed" > "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$SNa"_IGEref_out13.bed" -b "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
		# for TE-GENE fragments, the exon is the splice acceptor
		elif [ $a == "TE-GENE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_ACCEPTORS.bed" > "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$SNa"_IGEref_out13.bed" -b "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
		else
			c="N/A"
		fi
		if [[ $c = "" ]]; then c="."; fi
		# the loop aboce created the unix variable $c. here, this variable is appended to the main file.
		echo $line | awk -v c="$c" -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"c"\t"d}'
		rm "tmp."$SNa"_IGEref_out13.bed"
	done < $1 > $SNa"_IGEref_out15.tsv"
	echo " <-- done checking for breaks at splice sites at ... $(date)" >> $SNa"_PART5_"$logname".log"
}

split_TE_breakpoints()
{
	echo " --> start tyding up TE breakpoints at ... $(date)" >> $SNa"_PART5_"$logname".log"
	# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)
	cat $1 | while read line; do echo $line | awk '{print $8}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| sed 's/,$//'; done > $SNa"_newcolb"
	# for final output, get rid of original TE-breakpoint field
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13;}' $SNa"_IGEref_out15.tsv" > $SNa"_newcola"
	# append pooled TE breakpoints to final output
	paste $SNa"_newcola" $SNa"_newcolb" > $SNa"_chimericreads_final.tsv"
	rm $SNa"_newcola" $SNa"_newcolb"
	echo " <-- done tyding up TE breakpoints at ... $(date)" >> $SNa"_PART5_"$logname".log"
}

################################################################################
################################################################################
# execute:

# change to wd
cd $wd

# create .log file
logname=$(date | awk '{gsub(/\:/,"-",$5); print $4$3$2"_"$5}')
echo "====================" > $SNa"_PART5_"$logname".log"
echo "|| TEchim - PART5 || " >> $SNa"_PART5_"$logname".log"
echo "====================" >> $SNa"_PART5_"$logname".log"
echo "Parameters:" >> $SNa"_PART5_"$logname".log"
echo "Working directory:" "$wd" >> $SNa"_PART5_"$logname".log"
echo "Location of PART4 output:" "$path_to_PART4_output" >> $SNa"_PART5_"$logname".log"
echo "Sample name:" "$SNa" >> $SNa"_PART5_"$logname".log"
echo "Reference files:" "$REFpath$REFbase" >> $SNa"_PART5_"$logname".log"
echo "--------------------------------" >> $SNa"_PART5_"$logname".log"
echo " --> starting at ... $(date)" >> $SNa"_PART5_"$logname".log"
echo "--------------------------------" >> $SNa"_PART5_"$logname".log"

# create out10_combined.sorted.bed
cat $path_to_PART4_output"IGE_"$SNa*"/"$SNa*"_out10_breakpoints.bed" | bedtools sort -i - > $SNa"_IGEref_in10_combined.sorted.bed"

# add a gene tag to each read (i.e. genomic location). this is done using a gene-only version of the gtf file.
bedtools intersect -wa -a $SNa"_IGEref_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj -s > $SNa"_IGEref_out01_genetagged.tsv"

process_PART1_output $SNa"_IGEref_out01_genetagged.tsv"

combine_hits_of_each_TE $SNa"_IGEref_out03_TEbase.tsv"

check_for_splicesites $SNa"_IGEref_out12a_OUTPUT.tsv"

split_TE_breakpoints $SNa"_IGEref_out15.tsv"

echo " <-- all done at ... $(date)" >> $SNa"_PART5_"$logname".log"
