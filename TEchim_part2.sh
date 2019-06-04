#!/bin/bash

################################################################################
# TITLE: TEchim - PART 2
# VERSION: 0.2.0 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 04/06/2019 (dd/mm/yyyy)
# DESCRIPTION: This script combines separate samples and lanes, and generates
# summary .tsv files with all necessary data
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
path_to_PART1_output=/PATH/TO/PART1/OUTPUT/
SNa=NAME_OF_EXP
REFpath=/PATH/TO/REF/
REFbase=dmel625
################################################################################
################################################################################
# functions:

process_PART1_output()
{
	echo " --> start processing PART1 output at ... $(date)" >> $SNa"_PART2_"$logname".log"
	# separate | delimited field
	sed -e $'s/|/\t/g' $1 > $SNa"_out02_sepparated.tsv" && rm $SNa"_out01_genetagged.tsv"
	# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
	# and create column with BASE readname (no :A or :B)
		awk 'BEGIN {OFS = "\t"} {
		a = $5
		gsub(/_LTR/,"",a)
		c = substr($4, 1, length($4)-2)
		print $0"\t"a"\t"c
		}' < $SNa"_out02_sepparated.tsv" > $SNa"_out03pre_TEbase.tsv" && rm $SNa"_out02_sepparated.tsv"
	wc_predup=$(wc -l $SNa"_out03pre_TEbase.tsv" | awk '{print $1}')
	# remove duplicate reads (where both :A and :B version of the same reads were picked up
	awk '!seen[$21]++' $SNa"_out03pre_TEbase.tsv" > $SNa"_out03_TEbase.tsv" && rm $SNa"_out03pre_TEbase.tsv"
	echo " ---- $(wc -l $SNa"_out03_TEbase.tsv" | awk '{print $1}') reads ($wc_predup before removing duplicates)" >> $SNa"_PART2_"$logname".log"
	echo " <-- done processing PART1 output at ... $(date)" >> $SNa"_PART2_"$logname".log"
}

combine_hits_of_each_TE()
{
	echo " --> start combining reads at ... $(date)" >> $SNa"_PART2_"$logname".log"
	# create list of all TEs in dataset
	cut -f20 $1 | sort | uniq > $SNa"_out03a_uniqueTEs.tsv"
	# to create collection file, first make sure this file does not yet exist
	rm -f $SNa"_out12_OUTPUT.tsv"
	# loop through all TEs in the dataset
	while read TE
	do
		# grep TEs from combined data set
		awk -v t="$TE" '{if ($20 == t) print $0;}' $1 > "tmp."$SNa"_"$TE"_out04.tsv"
		if [ -s "tmp."$SNa"_"$TE"_out04.tsv" ]; then	
			###### SPLIT into TE-GENE and GENE-TE
			awk '{if ($7 == "TE-GENE") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv"
			# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
			# header: Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllBreakpoints(TE)
			if [ -s "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_out06_TE-GENE.tsv"
			fi
			awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv"
			if [ -s "tmp."$SNa"_"$TE"_out06_GENE-TE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_out06_GENE-TE.tsv"
			fi
			###### COMBINE TE-GENE and GENE-TE
			cat "tmp."$SNa"_"$TE"_out06_"*".tsv" > "tmp."$SNa"_"$TE"_out08.tsv"
			bedtools sort -i "tmp."$SNa"_"$TE"_out08.tsv" > "tmp."$SNa"_"$TE"_out09.tsv"
			# create list of all genes
			# first, replace "," with tab, then sort, then pick unique strings, then remove ""
			cut -f12 "tmp."$SNa"_"$TE"_out09.tsv" | tr ',' '\n' | sort | uniq | sed '/\./d' > "tmp."$SNa"_"$TE"_out11_genes.tsv"
			# this list of genes is used to add non-anchor hits to ANCHOR hits, if they overlap with the same gene and TE
			while read GENE
			do
				# create output file that contains for every gene all the insertions of that TE
				awk -v g="$GENE" '{if ($12 ~ g && g !~ $4) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$13"\t"$7"\t"$8"\t""mRNA";}' "tmp."$SNa"_"$TE"_out09.tsv" >> $SNa"_out12_OUTPUT.tsv"
			done < "tmp."$SNa"_"$TE"_out11_genes.tsv"
			rm -f "tmp."$SNa"_"$TE*
		fi
	done < $SNa"_out03a_uniqueTEs.tsv"
	# remove empty lines
	grep . $SNa"_out12_OUTPUT.tsv" > $SNa"_out12a_OUTPUT.tsv"
	rm $SNa"_out03a_uniqueTEs.tsv"
	rm $SNa"_out03_TEbase.tsv"
	echo " <-- done combining reads at ... $(date)" >> $SNa"_PART2_"$logname".log"
}

check_for_splicesites()
{
	echo " --> start checking for breaks at splice sites at ... $(date)" >> $SNa"_PART2_"$logname".log"
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
		d=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa -s | awk '{print $10}' | paste -sd ";" -)
		# check whether breakpoint on the genome overlaps with splice donor site (in the case of a GENE-TE fragment)
		# or with splice acceptor site (in the case of a TE-GENE fragment). the REF files _SPLICE_DONORS.bed and 
		# _SPLICE_ACCEPTORS.bed have a +/- 10nt, but the precise breakpoint is extracted here (which is the sequence
		# number of the last nucleotide that is still part of the exon for SPLICE_DONORS and the first exonic nt 
		# for SPLICE_ACCEPTORS.
		# for GENE-TE fragments, the exon is the splice donor
		if [ $a == "GENE-TE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_DONORS.bed" > "tmp."$SNa"_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b "tmp."$SNa"_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$SNa"_out14_ref_for_overlap.bed"
		# for TE-GENE fragments, the exon is the splice acceptor
		elif [ $a == "TE-GENE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_ACCEPTORS.bed" > "tmp."$SNa"_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b "tmp."$SNa"_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$SNa"_out14_ref_for_overlap.bed"
		else
			c="N/A"
		fi
		# the loop aboce created the unix variable $c. here, this variable is appended to the main file.
		echo $line | awk -v c="$c" -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"c"\t"d}'
		rm "tmp."$SNa"_out13.bed"
	done < $1 > $SNa"_out15.tsv"
	echo " <-- done checking for breaks at splice sites at ... $(date)" >> $SNa"_PART2_"$logname".log"
}

split_TE_breakpoints()
{
	echo " --> start tyding up TE breakpoints at ... $(date)" >> $SNa"_PART2_"$logname".log"
	# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)
	cat $1 | while read line; do echo $line | awk '{print $8}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| sed 's/,$//'; done > $SNa"_newcolb"
	# for final output, get rid of original TE-breakpoint field
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13;}' $SNa"_out15.tsv" > $SNa"_newcola"
	# append pooled TE breakpoints to final output
	paste $SNa"_newcola" $SNa"_newcolb" > $SNa"_chimericreads_final.tsv"
	rm $SNa"_newcola" $SNa"_newcolb"
	echo " <-- done tyding up TE breakpoints at ... $(date)" >> $SNa"_PART2_"$logname".log"
}

add_expression_levels()
{
	echo " --> start adding expression levels at ... $(date)" >> $SNa"_PART2_"$logname".log"
	list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}')
	while read line
	do
		# check if line in final .tsv has breakpoint near splice site ("." means no)
		if [[ $(echo $line | awk '{print $11}') == "." ]]; then
			# if breakpoint is not near splice site, copy the line
			echo $line
		else
			# get gene name
			input_gene=$(echo $line | awk '{print $1}')
			# get exons of input gene 				  | convert to bed 								  | sort			   | remove duplicate lines			> store as temporary file for this gene
			grep $input_gene $REFpath$REFbase"_EXONS.gtf" | awk '{print $1"\t"$4-1"\t"$5"\texon\t.\t"$7}' | bedtools sort -i - | awk '!seen[$2$3]++ {print $0}' > "tmp."$SNa"_"$input_gene".all_exons.bed"
			# create bedfile with chromosomal location of exon-intron junction
			echo $line | awk '{print $2"\t"$11-1"\t"$11"\t"$1"\t.\t"$3}' > "tmp."$SNa"_out1_breakpoint.bed"
			# create list with genomic distances of exons to chromosomal breakpoint -strand specific -print distance to (a) -show top 100
			bedtools closest -s -D a -k 100 -a "tmp."$SNa"_out1_breakpoint.bed" -b "tmp."$SNa"_"$input_gene".all_exons.bed" > "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv"
			# determine genomic region of exon that constitutes the chromosomal breakpoint (this assumes that the correct exon comes out on top after running bedtools closest) | convert to "region" syntax in samtools
			region1=$(head -n1 "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{print $7":"$8"-"$9}')
			# check if $region1 is empty
			if [[ -z $region1 ]]; then
				echo $line
			else
				# make sure collecting file does not yet exist
				rm -f "tmp."$SNa"_out5_TEreads"
				rm -f "tmp."$SNa"_out6_genereads"
				for snum in $list_of_snum
				do
					# make sure collecting file does not yet exist
					rm -f "tmp."$SNa"_"$snum"_out3_TEreads"
					rm -f "tmp."$SNa"_"$snum"_out4_genereads"
					list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2"\t"$3}' | grep $snum | awk '{print $2}')
					for lnum in $list_of_lanes
					do
						# check if STAR output bam has already been indexed, index if not
						if [[ -f $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam.bai" ]] ; then : ; else samtools index $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" ; fi
						# get number of reads inside region 1 where mate maps onto transposon
						samtools view $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v te="$(echo $line | awk '{print $5}')" '{ if($7 ~ te) {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$snum"_out3_TEreads"
						length_mode=$(samtools view $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" | sort | uniq -c | sort -rn | awk '{print $2}')
						# for the next step, both the strandedness of the gene and the type of fragment (GENE-TE or TE-GENE) determines the necessary steps
						if [[ $(echo $line | awk '{print $3}') == "+" ]]; then
							if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]; then
								# determine the nearest intron-exon junction beyond the breakpoint. if the gene is on the positive strand and the fragment is GENE-TE, then the next non-TE exon is located downstream (the smallest positive value in the bedtools closest output)
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13>0) {print$8}}' | head -n1)
								# first, get all reads that map in region1																											  | next, extract (and count) reads where mate maps beyond the next exon junction.
								samtools view $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$length_mode" '{ if($7 == "=" && $8>=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$snum"_out4_genereads"
							elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]; then
								# if fragment is TE-GENE, then the "next" intron-exon junction is upstream of the breakpoint ( the smallest NEGATIVE value in the bedtools closest output)
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13<0) {print$9}}' | head -n1)
								# first, get all reads that map in region1																	| next, extract (and count) reads where mate maps before the next exon junction.
								samtools view $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$length_mode" '{ if($7 == "=" && $8+lm<=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$snum"_out4_genereads"
							fi
						elif [[ $(echo $line | awk '{print $3}') == "-" ]]; then 
							if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]; then
								# if the gene is on the negative strand	and the fragment is GENE-TE, then the next non-TE exon is upstream of breakpoint (BUT: TAKE smallest POSITIVE value in bedtools closest output)
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13>0) {print$9}}' | head -n1)
								# first, get all reads that map in region1																	| next, extract (and count) reads where mate maps upstream to the next exon junction.
								samtools view $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$length_mode" '{ if($7 == "=" && ne!= "" && $8+lm<=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$snum"_out4_genereads"
							elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]; then
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13<0) {print$8}}' | head -n1)
								samtools view $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$length_mode" '{ if($7 == "=" && ne!= "" && $8>=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$snum"_out4_genereads"	
							fi
						fi
					done
					awk '{s+=$1} END {print s}' "tmp."$SNa"_"$snum"_out3_TEreads" >> "tmp."$SNa"_out5_TEreads"
					awk '{s+=$1} END {print s}' "tmp."$SNa"_"$snum"_out4_genereads" >> "tmp."$SNa"_out6_genereads"
				done			
				col1=$(paste -sd "|" "tmp."$SNa"_out5_TEreads")
				col2=$(paste -sd "|" "tmp."$SNa"_out6_genereads")
				echo $line | awk -v c1="$col1" -v c2="$col2" '{print $0"\tTE:"c1"\tGene:"c2}'
			fi
			rm -f "tmp."$SNa*
		fi
	done < $1 > $SNa"_chimericreads_final_withGENEreads.tsv"
	echo " <-- done adding expression levels at ... $(date)" >> $SNa"_PART2_"$logname".log"
}

################################################################################
################################################################################
# execute:

# change to wd
cd $wd

# create .log file
logname=$(date | awk '{gsub(/\:/,"-",$5); print $4$3$2"_"$5}')
echo "====================" > $SNa"_PART2_"$logname".log"
echo "|| TEchim - PART2 || " >> $SNa"_PART2_"$logname".log"
echo "====================" >> $SNa"_PART2_"$logname".log"
echo "Parameters:" >> $SNa"_PART2_"$logname".log"
echo "Working directory:" "$wd" >> $SNa"_PART2_"$logname".log"
echo "Location of PART1 output:" "$path_to_PART1_output" >> $SNa"_PART2_"$logname".log"
echo "Sample name:" "$SNa" >> $SNa"_PART2_"$logname".log"
echo "Reference files:" "$REFpath$REFbase" >> $SNa"_PART2_"$logname".log"
echo "--------------------------------" >> $SNa"_PART2_"$logname".log"
echo " --> starting at ... $(date)" >> $SNa"_PART2_"$logname".log"
echo "--------------------------------" >> $SNa"_PART2_"$logname".log"

# create out10_combined.sorted.bed
cat $path_to_PART1_output$SNa*"/"$SNa*"_out10_breakpoints.bed" | bedtools sort -i - > $SNa"_in10_combined.sorted.bed"

# add a gene tag to each read (i.e. genomic location). this is done using a gene-only version of the gtf file.
bedtools intersect -wa -a $SNa"_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj -s > $SNa"_out01_genetagged.tsv"

process_PART1_output $SNa"_out01_genetagged.tsv"

combine_hits_of_each_TE $SNa"_out03_TEbase.tsv"

check_for_splicesites $SNa"_out12a_OUTPUT.tsv"

split_TE_breakpoints $SNa"_out15.tsv"

# optional. Only recommended AFTER filtering.
#add_expression_levels $SNa"_chimericreads_final.tsv"

echo " <-- all done at ... $(date)" >> $SNa"_PART2_"$logname".log"
