#!/bin/bash

################################################################################
# TITLE: TEchim - PART 2-5
# VERSION: 0.2.1 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 05/06/2019 (dd/mm/yyyy)
# DESCRIPTION: This script combines parts 2-5. 
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
wd=~							# working directory (no trailing "/")
path_to_PART1_output=/PATH/TO/OUTPUT/	# with trailing "/"
SNa=NAME_OF_EXP					# sample name	
nc=1							# number of cores
REFpath=/PATH/TO/REF/			# same path as in _buildREF (with trailing "/")
REFbase=dmel625					# same name as in _buildREF
################################################################################
################################################################################
# functions:

process_P1out_TE()
{
	cd $wd
	echo " --> start processing PART1 output for TE at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
	cat $1 | bedtools sort -i - > $SNa"_in10_combined.sorted.bed"
	bedtools intersect -wa -a $SNa"_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj -s > $SNa"_out01_genetagged.tsv"
	# separate | delimited field
	sed -e $'s/|/\t/g' $SNa"_out01_genetagged.tsv" > $SNa"_out02_sepparated.tsv" && rm $SNa"_out01_genetagged.tsv"
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
	echo " ---- $(wc -l $SNa"_out03_TEbase.tsv" | awk '{print $1}') TE reads ($wc_predup before removing duplicates)" >> $SNa"_PART2to5_"$logname".log"
	echo " <-- done processing PART1 output for TE at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

combine_hits_of_each_TE()
{
	cd $wd
	echo " --> start combining TE reads at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
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
				awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$13"\t"$7"\t"$8"\t""mRNA";}' "tmp."$SNa"_"$TE"_out09.tsv" >> $SNa"_out12_OUTPUT.tsv"
			done < "tmp."$SNa"_"$TE"_out11_genes.tsv"
			rm -f "tmp."$SNa"_"$TE*
		fi
	done < $SNa"_out03a_uniqueTEs.tsv"
	# remove empty lines
	grep . $SNa"_out12_OUTPUT.tsv" > $SNa"_out12a_OUTPUT.tsv"
	rm $SNa"_out03a_uniqueTEs.tsv"
	rm $SNa"_out03_TEbase.tsv"
	echo " <-- done combining TE reads at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

check_for_TE_splicesites()
{
	cd $wd
	echo " --> start checking for TE breaks at splice sites at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
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
		if [[ $c = "" ]]; then c="."; fi
		# the loop aboce created the unix variable $c. here, this variable is appended to the main file.
		echo $line | awk -v c="$c" -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"c"\t"d}'
		rm "tmp."$SNa"_out13.bed"
	done < $1 > $SNa"_out15.tsv"
	echo " <-- done checking for TE breaks at splice sites at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

split_TE_breakpoints()
{
	cd $wd
	echo " --> start tyding up TE breakpoints at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
	# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)
	cat $1 | while read line; do echo $line | awk '{print $8}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| sed 's/,$//'; done > $SNa"_newcolb"
	# for final output, get rid of original TE-breakpoint field
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13;}' $SNa"_out15.tsv" > $SNa"_newcola"
	# append pooled TE breakpoints to final output
	paste $SNa"_newcola" $SNa"_newcolb" > $SNa"_TE_chimericreads_final.tsv"
	rm $SNa"_newcola" $SNa"_newcolb"
	echo " <-- done tyding up TE breakpoints at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

add_expression_levels()
{
	cd $wd
	echo " --> start adding expression levels at ... $(date)" >> $SNa"_PART2_"$logname".log"
	# get length of aligned files - here only the first 1000 reads from the first alignment file is consiedered
	# get the name of the first sample
	first_file_name=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | head -n1)
	# get first 1000 aligned reads																												| print length of sequences		 | sort in descending order and | take top
	length_mode=$(samtools view $path_to_PART1_output$first_file_name"/"$first_file_name"_STAR/"$first_file_name"_out4_Aligned.sortedByCoord.out.bam" | head -n1000 | awk '{a=length($10); print a}' | sort -rn | head -n1)
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
	done < $1 > $SNa"_TE_chimericreads_final_withGENEreads.tsv"
	echo " <-- done adding expression levels at ... $(date)" >> $SNa"_PART2_"$logname".log"
}

find_matching_IGEs()
{
	cd $wd
	echo " --> start finding matching IGEs at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
	# randomly subsample filtered CDS (here: take 10%)
	tenpc_CdS=$(wc -l $REFpath$REFbase".CDS_for_IGE.bed" | awk '{OFMT = "%.0f"; print $1/10}')
	sort -R $REFpath$REFbase".CDS_for_IGE.bed" | head -n $tenpc_CdS > $SNa"_inputGENEs.bed"
	# assess how many Samples and Lanes exist - run bedtools mutlicov on each sample/lane
	list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}')
	for snum in $list_of_snum
	do
		list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2"\t"$3}' | grep $snum | awk '{print $2}')
		for lnum in $list_of_lanes
		do
			# check if .bam file is indexed, if not then index it
			if [[ ! -f $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam.bai" ]]; then samtools index $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam"; fi
			# quantify coverage of TEs in each sample/lane - this line will print the number
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $REFpath"TEs.fa.bed" -s | awk '{print $7}' > "tmp."$SNa"_"$snum"_"$lnum"_out31_TEcoverage_perSample.tsv"
			# this line will print the TE name - only has to be done once (hence the if statement)
			if [[ ! -f "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv" ]]; then bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $REFpath"TEs.fa.bed" -s | awk '{print $4}' > "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv"; fi
			# quantify coverage of CDS (subsample) in each sample/lane - this line will print the number
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $SNa"_inputGENEs.bed" -s | awk '{print $7}' > "tmp."$SNa"_"$snum"_"$lnum"_out31_GENEcoverage_perSample.tsv"
			# this line will print the GENE names - only has to be done once (hence the if statement)
			if [[ ! -f "tmp."$SNa"_out31a_GENEcoverage_GENEnames.tsv" ]]; then bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $SNa"_inputGENEs.bed" -s | awk '{print $4}' > "tmp."$SNa"_out31a_GENEcoverage_GENEnames.tsv"; fi	
		done
	done
	# combine all TE coverage numbers					  | calculate average, round to full integer
	paste "tmp."$SNa"_"*"_out31_TEcoverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$SNa"_out31b_TEcoverage_averages.tsv"
	# combine averages with TE names
	paste "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv" "tmp."$SNa"_out31b_TEcoverage_averages.tsv" | sort -k2n > "tmp."$SNa"_out32_TEcoverage_averages.tsv"
	# combine all GENE coverage numbers						| calculate average, round to full integer
	paste "tmp."$SNa"_"*"_out31_GENEcoverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$SNa"_out31b_GENEcoverage_averages.tsv"
	# combine averages with GENE names
	paste "tmp."$SNa"_out31a_GENEcoverage_GENEnames.tsv" "tmp."$SNa"_out31b_GENEcoverage_averages.tsv" | sort -k2n > "tmp."$SNa"_out32_GENEcoverage_averages.tsv"
	# loop through all TEs
	while read TEline
	do
		coverage=$(echo $TEline | awk '{print $2}')
		# find the gene that matches the expression level of each TE									  | print matching TE-GENE pairs
		awk -v v1="$coverage" '$2>v1 {print $1"\t"$2; exit}' "tmp."$SNa"_out32_GENEcoverage_averages.tsv" | awk -v te="$TEline" '{print te"\t"$0}'
	done < "tmp."$SNa"_out32_TEcoverage_averages.tsv" > "tmp."$SNa"_out35_TEmatchedCDS.tsv"
	# left-join TE-GENE pairs with original GENE bed (this is to get the strand of the according CDS's)																				 | print a bed file with each UNIQUE CDS (strand-specific!)
	join <(awk '{print $3"\t"$1"\t"$2"\t"$4}' "tmp."$SNa"_out35_TEmatchedCDS.tsv" | sort -t $'\t') <(awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $SNa"_inputGENEs.bed" | sort -t $'\t') | awk '!seen[$1]++ {print $5"\t"$6"\t"$7"\t"$1"\t.\t"$8}' > $SNa"_IGEs.bed"
	# also print a lookup table to link CDS's to the matching TE
	join <(awk '{print $3"\t"$1"\t"$2"\t"$4}' "tmp."$SNa"_out35_TEmatchedCDS.tsv" | sort -t $'\t') <(awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $SNa"_inputGENEs.bed" | sort -t $'\t') | awk '{print $1"\t"$4"\t"$2"\t"$3}' > $SNa"_IGE_TE_lookup.tsv"
	rm "tmp."$SNa*
	echo " <-- done finding matching IGEs at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

create_IGE_reference()
{
	cd $wd
	echo " --> start creating IGE reference files at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
	mkdir $SNa"_IGEref_CONTAINER"
	cd $SNa"_IGEref_CONTAINER"
	# get fasta for CDS's
	bedtools getfasta -fi $REFpath$REFbase".fa" -bed $1 -s -name > $SNa"_IGEref.TEs.fa"
	# copy original reference genome (this is necessary because RepeatMasker changes current files)
	cp $REFpath$REFbase".fa" "./tmp."$REFbase".fa"
	# run RepeatMasker using newly created CDS fasta file as "TE", then remove unnecessary output
	RepeatMasker -lib $SNa"_IGEref.TEs.fa" -no_is -nolow -s -pa $nc "./tmp."$REFbase".fa"
	rm "./tmp."$REFbase".fa"
	rm $REFbase".fa.cat.gz"
	rm $REFbase".fa.ori.out"
	rm $REFbase".fa.out"
	rm $REFbase".fa.tbl"
	# output is a "clean" reference genome, now void of IGE sequences
	mv $REFbase".fa.masked" $SNa"_IGEref.clean.noTEs.fa"
	# create FASTA of IGEs where each IGE has chromosome name ">TEchr_..."
	awk '{if ($1 ~ ">") {gsub(/>/,""); print ">TEchr_"$1"\t"$2"\t"$3} else {print $0}}' $SNa"_IGEref.TEs.fa" > $SNa"_IGEref.clean.onlyTEs.fa"
	# combine IGE-cleaned reference genome with IGE-fasta file	
	cat $SNa"_IGEref.clean.noTEs.fa" $SNa"_IGEref.clean.onlyTEs.fa" > $SNa"_IGEref.clean.fa"
	# generate STAR genome index
	mkdir "STAR_"$SNa"_IGEref"
	STAR --runMode genomeGenerate --genomeFastaFiles $SNa"_IGEref.clean.fa" --genomeDir "STAR_"$SNa"_IGEref" --runThreadN $nc
	mv Log.out "STAR_"$SNa"_IGEref/."
	# generate BLAST databases
	makeblastdb -dbtype nucl -in $SNa"_IGEref.clean.onlyTEs.fa"
	makeblastdb -dbtype nucl -in $SNa"_IGEref.clean.noTEs.fa"
	cd $wd
	echo " <-- done creating IGE reference files at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

align_IGEref_and_filter()
{
	cd $wd
	if [ ! -d "IGE_"$SNa"_"$SNo"_"$LNo ]; then mkdir "IGE_"$SNa"_"$SNo"_"$LNo; fi
	cd "IGE_"$SNa"_"$SNo"_"$LNo
	echo "====================" > $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "|| TEchim - PART4 || " >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "====================" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " --> start mapping on IGE reference at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# run STAR in chimera-mode
	STAR --runThreadN $nc \
		--genomeDir $wd"/"$SNa"_IGEref_CONTAINER/STAR_"$SNa"_IGEref" \
		--readFilesIn $1 $2 \
		--chimSegmentMin 20 \
		--chimOutType WithinBAM \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix $SNa"_IGEref_"$SNo"_"$LNo$"_out4_"	
	# extract only hits that cross IGE-GENE breakpoints. the awk commands remove IGE-IGE reads
	samtools view $SNa"_IGEref_"$SNo"_"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | grep TEchr_ | awk '($7 != "=")' > $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam"
	mkdir $SNa"_IGEref_"$SNo"_"$LNo"_STAR"
	mv $SNa"_IGEref_"$SNo"_"$LNo$"_out4"* $SNa"_IGEref_"$SNo"_"$LNo"_STAR"/.
	if [ ! -s $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam" ]; then
		echo " #### ERROR: file does not have any TE-GENE brakpoint spanning reads!" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
		echo " --> exited at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
		exit
	fi
	echo " ------ sample contains $(wc -l $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam" | awk '{print $1}') reads that span gene-TE breakpoint." >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " <-- done with mapping on IGE reference at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "--------------------------------" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	cd $wd
}

blast_on_longreads ()
{
	cd $wd"/IGE_"$SNa"_"$SNo"_"$LNo
	echo " --> start BLAST at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# extract readnames and remove duplicates
	awk '{print $1}' $1 | sort | uniq > $SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt"
	# combine readnames with long sequences stored in the lookup file
	# LANG=en_EN is a bug-fix to make sort compatible with join
	LANG=en_EN join -1 1 -2 1 <(sort $SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt") <(zcat < $2 ) > $SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv" && rm $SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt"
	# the following while loop will add data to file using ">>". just as safety
	# precaution, this line makes sure no data with such a name exists
	rm -f $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
	# loop through long sequences. two blast searches are performed, one on a
	# reference genome without TE sequences and another one on all TE sequences.
	# only the top hit is reported.
	while IFS=$' ' read -r readname sequence
	do	
		echo "${readname}" > var1
		echo "${sequence}" > var2
		echo "${sequence}" | blastn -db $wd"/"$SNa"_IGEref_CONTAINER/"$SNa"_IGEref.clean.noTEs.fa" -outfmt "6 qstart qend sseqid sstart send sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var3		
		if ! [ -s var3 ]; then echo "no-hit" > var3; fi
		echo "${sequence}" | blastn -db $wd"/"$SNa"_IGEref_CONTAINER/"$SNa"_IGEref.clean.onlyTEs.fa" -outfmt "6 qstart qend sseqid sstart send slen sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var4
		if ! [ -s var4 ]; then echo "no-hit" > var4; fi
		paste var1 var2 var3 var4 >> $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
		rm var*
	done < $SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv"
	# remove reads that did not give BLAST result
	grep -v no-hit $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv" > $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv"
	echo " ------ BLAST results: $(wc -l $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" | awk '{print $1}') hits ($(grep no-hit $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv" | wc -l | awk '{print $1}') no hit)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm -f $SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv"
	rm -f $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
	echo " <-- done with BLAST at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "--------------------------------" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam"
	cd $wd
}

create_summary_table ()
{
	cd $wd"/IGE_"$SNa"_"$SNo"_"$LNo
	echo " --> start creating summary table at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# determine whether the section that maps to the genome is the 5' or the 3' end
	# of the mRNA section:
	# 5'-######|GENE|TE|######-3' => PART5
	# 5'-######|TE|GENE|######-3' => PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; if (a < b) {print "PART5"} else {print "PART3"}}' < $1 > tmpfile.genepart
	# determine the precise breakpoint on the chromosome. this depends on  whether
	# the chromosomal part is PART5 or PART3
	awk 'BEGIN {OFS = "\t"} {a = $5; {print a}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.chr
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d-1} else {print c-1}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.start
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d} else {print c}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.end
	# determine the precise breakpoint on the TE. this depends on  whether the TE
	# part is PART5 or PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $13 ; d = $12 ; if (a < b) {print d} else {print c}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.TE
	# determine the overlap between the two mapped sections of the long read
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $4 ; d = $10 ; if (a < b) {print b-c-1} else {print a-d-1}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.uncertainty
	if [[ $stranded = "0" ]]; then
		awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
			a = $11
			gsub(/TEchr_/,"",a)
			if ($8 == "plus") {if ($15=="plus") {b = "forward"} else {b = "reverse"}} else {if ($15 == "plus") { b = "reverse" } else { b = "forward" }}
			if ($3 < $9) {print $1"|"a"|"b"|GENE-TE|S"s"|L"l} else {print $1"|"a"|"b"|TE-GENE|S"s"|L"l}
			}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.readname
			awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "+"}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.strand
	else
	awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
		a = $11
		gsub(/TEchr_/,"",a)
		b = $15
		c = $3
		d = $9
		e = $1
		if (c < d) {print e"|"a"|"b"|GENE-TE|S"s"|L"l} else {print e"|"a"|"b"|TE-GENE|S"s"|L"l}
		}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.readname
		awk 'BEGIN {OFS = "\t"} {a = $8 ; if (a == "plus") {print "+"} else {print "-"}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.strand
	fi
	awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "."}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.score
	paste -d'|' tmpfile.readname tmpfile.breakpoint.TE tmpfile.uncertainty > tmpfile.readname.extended
	paste -d'\t' tmpfile.chr tmpfile.breakpoint.chr.start tmpfile.breakpoint.chr.end tmpfile.readname.extended tmpfile.score tmpfile.breakpoint.chr.strand > $SNa"_IGEref_"$SNo"_"$LNo$"_out10_breakpoints.bed"
	#paste -d'\t' $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" tmpfile.breakpoint.chr.end tmpfile.breakpoint.TE tmpfile.uncertainty > $SNa"_IGEref_"$SNo"_"$LNo$"_out11_combined_results.tsv" &&
	rm tmpfile.*
	echo " <-- all done at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "================================" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv"
	cd $wd
}

process_P1out_IGE()
{
	cd $wd
	echo " --> start processing PART1 output for IGE at ... $(date)" >> $SNa"_PART2to5_"$logname".log"	
	cat $1 | bedtools sort -i - > $SNa"_IGEref_in10_combined.sorted.bed"
	bedtools intersect -wa -a $SNa"_IGEref_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj -s > $SNa"_IGEref_out01_genetagged.tsv"
	# filter reads inside same gene
	awk '{if ($4 !~ $10) print $0}' $SNa"_IGEref_out01_genetagged.tsv" > $SNa"_IGEref_out01a_filtered_genetagged.tsv"
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
	echo " ---- $(wc -l $SNa"_IGEref_out03_TEbase.tsv" | awk '{print $1}') reads ($wc_predup before removing duplicates)" >> $SNa"_PART2to5_"$logname".log"
	echo " <-- done processing PART1 output for IGE at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

combine_hits_of_each_IGE()
{
	cd $wd
	echo " --> start combining IGE reads at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
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
			# SPLIT into TE-GENE and GENE-TE
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
			# COMBINE TE-GENE and GENE-TE
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
	echo " <-- done combining IGE reads at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

check_for_IGE_splicesites()
{
	echo " --> start checking for IGE breaks at splice sites at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
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
	echo " <-- done checking for IGE breaks at splice sites at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}

split_IGE_breakpoints()
{
	cd $wd
	echo " --> start tyding up IGE breakpoints at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
	# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)
	cat $1 | while read line; do echo $line | awk '{print $8}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| sed 's/,$//'; done > $SNa"_IGE_newcolb"
	# for final output, get rid of original TE-breakpoint field
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13;}' $SNa"_IGEref_out15.tsv" > $SNa"_IGE_newcola"
	# append pooled TE breakpoints to final output
	paste $SNa"_IGE_newcola" $SNa"_IGE_newcolb" > $SNa"_IGEref_chimericreads_final.tsv"
	rm $SNa"_IGE_newcola" $SNa"_IGE_newcolb"
	echo " <-- done tyding up IGE breakpoints at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
}


################################################################################
################################################################################
# execute:

# change to wd
cd $wd
# create .log file
logname=$(date | awk '{gsub(/\:/,"-",$5); print $4$3$2"_"$5}')
echo "======================" > $SNa"_PART2to5_"$logname".log"
echo "|| TEchim - PART2-5 || " >> $SNa"_PART2to5_"$logname".log"
echo "======================" >> $SNa"_PART2to5_"$logname".log"
echo "Parameters:" >> $SNa"_PART2to5_"$logname".log"
echo "Working directory:" "$wd" >> $SNa"_PART2to5_"$logname".log"
echo "Location of PART1 output:" "$path_to_PART1_output" >> $SNa"_PART2to5_"$logname".log"
echo "Sample name:" "$SNa" >> $SNa"_PART2to5_"$logname".log"
echo "Reference files:" "$REFpath$REFbase" >> $SNa"_PART2to5_"$logname".log"
echo "--------------------------------" >> $SNa"_PART2to5_"$logname".log"
echo " --> starting at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
echo "--------------------------------" >> $SNa"_PART2to5_"$logname".log"

# PART2:
(process_P1out_TE $path_to_PART1_output$SNa*"/"$SNa*"_out10_breakpoints.bed"
combine_hits_of_each_TE $SNa"_out03_TEbase.tsv"
check_for_TE_splicesites $SNa"_out12a_OUTPUT.tsv"
split_TE_breakpoints $SNa"_out15.tsv"
# optional. Only recommended AFTER filtering.
#add_expression_levels $SNa"_TE_chimericreads_final.tsv"
) &

# PART3:
find_matching_IGEs
create_IGE_reference $SNa"_IGEs.bed"

# PART4:
list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}')
for SNo in $list_of_snum
do
	list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2"\t"$3}' | grep $snum | awk '{print $2}')
	for LNo in $list_of_lanes
	do
		(align_IGEref_and_filter $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_out3_1.fasta" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_out3_2.fasta" 
		blast_on_longreads $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_LOOKUP.sorted.tsv.gz"
		create_summary_table $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv") &
	done
done
wait

# PART5:
process_P1out_IGE $wd"/IGE_"$SNa*"/"$SNa*"_out10_breakpoints.bed"
combine_hits_of_each_IGE $SNa"_IGEref_out03_TEbase.tsv"
check_for_IGE_splicesites $SNa"_IGEref_out12a_OUTPUT.tsv"
split_IGE_breakpoints $SNa"_IGEref_out15.tsv"
echo " <-- all done at ... $(date)" >> $SNa"_PART2to5_"$logname".log"
