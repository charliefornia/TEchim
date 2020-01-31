#!/bin/bash

################################################################################
# TITLE: TEchim - PART 2
# VERSION: 0.2.3 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 24/06/2019 (dd/mm/yyyy)
# DESCRIPTION: This script combines output from PART1
################################################################################

################################################################################
# REQUIREMENTS:
# - STAR (https://github.com/alexdobin/STAR)
# - samtools
# - bedtools
# - blast
################################################################################

################################################################################
################################################################################
# set parameters
wd=$(pwd)						# working directory (no trailing "/")
path_to_PART1_output=$wd"/"		# with trailing "/"
SNa=
nc=1							# number of cores (default:1)
REFpath=						# if left empty then REFpath from PART1 is used
################################################################################
################################################################################
# functions:

get_variables()
{
	# assign REFpath parameter
	if [ -z $REFpath ]; then
		if [ -e $path_to_PART1_output"."$SNa"_refpath" ]; then
			REFpath=$(cat $path_to_PART1_output"."$SNa"_refpath")
		else
			echo " #### ERROR: path to output from PART1 is corrupt - no file named .""$SNa""_refpath"
			exit
		fi
	fi
	# assign REFbase parameter
	if [ -e $REFpath"REFERENCE_basename" ]; then
		REFbase=$(cat $REFpath"REFERENCE_basename")
	else
		echo " #### ERROR: reference path is corrupt - no file named REFERENCE_basename"
		exit
	fi
	# check strandedness of input FASTQ
	if [ -e $path_to_PART1_output"."$SNa"_strandedness" ]; then
		stranded=$(cat $path_to_PART1_output"."$SNa"_strandedness")
	else
		echo " #### ERROR: path to output from PART1 is corrupt - no file named .""$SNa""_strandedness"
		exit
	fi
	# check readlength of input FASTQ
	if [ -e $path_to_PART1_output"."$SNa"_fastalength" ]; then
		fastalength=$(cat $path_to_PART1_output"."$SNa"_fastalength")
	else
		echo " #### ERROR: path to output from PART1 is corrupt - no file named .""$SNa""_fastalength"
		exit
	fi
}

write_logfile()
{
	logname=$(date | awk '{gsub(/\:/,"-",$5); print $4$3$2"_"$5}')
	echo "======================" > $wd"/"$SNa"_PART2_"$logname".log"
	echo "|| TEchim - PART2-5 || " >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "======================" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "Parameters:" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "Working directory:" "$wd" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "Location of PART1 output:" "$path_to_PART1_output" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "Sample name:" "$SNa" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "Reference files:" "$REFpath$REFbase" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "--------------------------------" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo " --> starting at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "--------------------------------" >> $wd"/"$SNa"_PART2_"$logname".log"
}

process_P1out_TE()
{
	cd $wd
	echo " --> start processing PART1 output for TE at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
	cat $path_to_PART1_output$SNa*"/"$SNa*"_out10_breakpoints.bed" | bedtools sort -i - > $SNa"_in10_combined.sorted.bed"
	if [ $stranded = "0" ]; then
		bedtools intersect -wa -a $SNa"_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj > $SNa"_out01_genetagged.tsv"
	else
		bedtools intersect -wa -a $SNa"_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj -s > $SNa"_out01_genetagged.tsv"
	fi
	
	# Also extract those reads that were not inside annotated gene
	#grep -Fvf <(tr "|" "\t" < $SNa"_in10_combined.sorted.bed" | awk '{print $4}') <(cat $SNa"_out01_genetagged.tsv") > $SNa"_out01b_ReadsOutsideGenes.tsv"
	
	# separate | delimited field
	tr '|' '\t' < $SNa"_out01_genetagged.tsv" > $SNa"_out02_sepparated.tsv" && rm $SNa"_out01_genetagged.tsv"
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
	rm $SNa"_in10_combined.sorted.bed"
	echo " ---- $(wc -l $SNa"_out03_TEbase.tsv" | awk '{print $1}') TE reads ($wc_predup before removing duplicates)" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo " <-- done processing PART1 output for TE at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
}

combine_hits_of_each_TE()
{
	cd $wd
	echo " --> start combining TE reads at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
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
			# header: Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllSamples|AllBreakpoints(TE)
			if [ -s "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,6,7,17,8,10 -o distinct,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse,collapse -d 20 > "tmp."$SNa"_"$TE"_out06_TE-GENE.tsv"
			fi
			awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv"
			if [ -s "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,6,7,17,8,10 -o distinct,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse,collapse -d 20 > "tmp."$SNa"_"$TE"_out06_GENE-TE.tsv"
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
				if [ $stranded = "0" ]; then
					awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$7"\t"$8"\tunstranded\t"$13"\t"$14;}' "tmp."$SNa"_"$TE"_out09.tsv" >> $SNa"_out12_OUTPUT.tsv"
				else
					awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$7"\t"$8"\tmRNA\t"$13"\t"$14;}' "tmp."$SNa"_"$TE"_out09.tsv" >> $SNa"_out12_OUTPUT.tsv"
				fi
			done < "tmp."$SNa"_"$TE"_out11_genes.tsv"
			rm -f "tmp."$SNa"_"$TE*
		fi
	done < $SNa"_out03a_uniqueTEs.tsv"
	# remove empty lines
	grep . $SNa"_out12_OUTPUT.tsv" > $SNa"_out12a_OUTPUT.tsv"
	rm $SNa"_out12_OUTPUT.tsv"
	rm $SNa"_out03a_uniqueTEs.tsv"
	rm $SNa"_out03_TEbase.tsv"
	echo " <-- done combining TE reads at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
}

check_for_TE_splicesites()
{
	cd $wd
	echo " --> start checking for TE breaks at splice sites at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
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
		if [ $stranded = "0" ]; then
			d=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa | awk '{print $10}' | paste -sd ";" -)
		else
			d=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa -s | awk '{print $10}' | paste -sd ";" -)
		fi
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
		if [[ $c = "" ]]; then c="."; fi
		# the loop aboce created the unix variable $c. here, this variable is appended to the main file.
		echo $line | awk -v c="$c" -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"c"\t"d}'
		rm "tmp."$SNa"_out13.bed"
	done < $1 > $SNa"_out15.tsv"
	rm $1
	echo " <-- done checking for TE breaks at splice sites at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
}

split_TE_breakpoints()
{
	cd $wd
	echo " --> start tyding up TE breakpoints at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
	# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)	
	cat $1 | while read line; do echo $line | awk '{print $12}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| awk '{print $0"\n"}' | awk 'NF'; done > $SNa"_newcolb"
	cat $1 | while read line; do echo $line | awk '{print $11}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| awk '{print $0"\n"}' | awk 'NF'; done > $SNa"_newcolc"
	# for final output, get rid of original TE-breakpoint field
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$13,$14}' $SNa"_out15.tsv" > $SNa"_newcola"
	# append pooled TE breakpoints to final output
	paste $SNa"_newcola" $SNa"_newcolb" $SNa"_newcolc" > $SNa"_TE_chimericreads_final.tsv"
	rm $SNa"_newcola" $SNa"_newcolb" $SNa"_newcolc"
	rm $1
	echo " <-- done tyding up TE breakpoints at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
}

add_expression_levels()
{
	cd $wd
	echo " --> start adding expression levels at ... $(date)" >> $SNa"_PART2_"$logname".log"
	list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
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
				for SNo in $list_of_SNo
				do
					# make sure collecting file does not yet exist
					rm -f "tmp."$SNa"_"$SNo"_out3_TEreads"
					rm -f "tmp."$SNa"_"$SNo"_out4_genereads"
					list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
					for LNo in $list_of_LNo
					do
						# check if STAR output bam has already been indexed, index if not
						if [[ -f $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam.bai" ]] ; then : ; else samtools index $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" ; fi
						# get number of reads inside region 1 where mate maps onto transposon
						samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v te="$(echo $line | awk '{print $5}')" '{ if($7 ~ te) {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$SNo"_out3_TEreads"
						# for the next step, both the strandedness of the gene and the type of fragment (GENE-TE or TE-GENE) determines the necessary steps
						if [[ $(echo $line | awk '{print $3}') == "+" ]]; then
							if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]; then
								# determine the nearest intron-exon junction beyond the breakpoint. if the gene is on the positive strand and the fragment is GENE-TE, then the next non-TE exon is located downstream (the smallest positive value in the bedtools closest output)
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13>0) {print$8}}' | head -n1)
								# first, get all reads that map in region1																											  | next, extract (and count) reads where mate maps beyond the next exon junction.
								samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$fastalength" '{ if($7 == "=" && $8>=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$SNo"_out4_genereads"
							elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]; then
								# if fragment is TE-GENE, then the "next" intron-exon junction is upstream of the breakpoint ( the smallest NEGATIVE value in the bedtools closest output)
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13<0) {print$9}}' | head -n1)
								# first, get all reads that map in region1																	| next, extract (and count) reads where mate maps before the next exon junction.
								samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$fastalength" '{ if($7 == "=" && $8+lm<=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$SNo"_out4_genereads"
							fi
						elif [[ $(echo $line | awk '{print $3}') == "-" ]]; then 
							if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]; then
								# if the gene is on the negative strand	and the fragment is GENE-TE, then the next non-TE exon is upstream of breakpoint (BUT: TAKE smallest POSITIVE value in bedtools closest output)
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13>0) {print$9}}' | head -n1)
								# first, get all reads that map in region1																	| next, extract (and count) reads where mate maps upstream to the next exon junction.
								samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$fastalength" '{ if($7 == "=" && ne!= "" && $8+lm<=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$SNo"_out4_genereads"
							elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]; then
								next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13<0) {print$8}}' | head -n1)
								samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" -v lm="$fastalength" '{ if($7 == "=" && ne!= "" && $8>=ne && $6 == lm"M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_"$SNo"_out4_genereads"	
							fi
						fi
					done
					awk '{s+=$1} END {print s}' "tmp."$SNa"_"$SNo"_out3_TEreads" >> "tmp."$SNa"_out5_TEreads"
					awk '{s+=$1} END {print s}' "tmp."$SNa"_"$SNo"_out4_genereads" >> "tmp."$SNa"_out6_genereads"
				done			
				col1=$(paste -sd "|" "tmp."$SNa"_out5_TEreads")
				col2=$(paste -sd "|" "tmp."$SNa"_out6_genereads")
				echo $line | awk -v c1="$col1" -v c2="$col2" '{print $0"\tTE:"c1"\tGene:"c2}'
			fi
			rm -f "tmp."$SNa*
		fi
	done < $1 > $SNa"_TE_chimericreads_final_withGENEreads.tsv"
	echo " <-- done adding expression levels at ... $(date)" >> $SNa"_PART2_"$logname".log"
}



################################################################################
################################################################################
# execute:

# change to wd
cd $wd

get_variables
write_logfile
process_P1out_TE
combine_hits_of_each_TE $wd"/"$SNa"_out03_TEbase.tsv"
if [ $stranded != "0" ]; then
	check_for_TE_splicesites $wd"/"$SNa"_out12a_OUTPUT.tsv"
	split_TE_breakpoints $wd"/"$SNa"_out15.tsv"
	# optional. Only recommended AFTER filtering.
	#add_expression_levels $SNa"_TE_chimericreads_final.tsv"
else
	split_TE_breakpoints $wd"/"$SNa"_out12a_OUTPUT.tsv"
fi
echo " <-- all done at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
