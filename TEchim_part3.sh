#!/bin/bash

################################################################################
# TITLE: TEchim - PART3 - detection of TE-gene chimera in RNA-seq data
# VERSION: 0.1.2 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 17/05/2019 (dd/mm/yyyy)
# DESCRIPTION:
################################################################################

################################################################################
# REQUIREMENTS:
# - samtools
# - bedtools
################################################################################

################################################################################
# set parameters
wd=~/Dropbox/CloudDesktop/TEchim_cloud/ANALYSIS/
path_to_TEchimPART2_output=~/Documents/2019MAY_TEscoex_longRNA/fromHPC/
PART2_output=~/Dropbox/CloudDesktop/TEchim_cloud/2018MARCH_TEchim_chimericreads_final.tsv
SNa=2018MARCH_TEchim
NumSam=6
NumLan=2
REF=~/Dropbox/CloudDesktop/REF_cloud/
REFbase=dmel625_v04
################################################################################

# NOTE: this is assuming that read length is 60nt

# change to wd
cd $wd
while read line
do
	# check if line in final .tsv has breakpoint near splice site ("." means no)
	if [[ $(echo $line | awk '{print $11}') == "." ]]
	then
	# if breakpoint is not near splice site, copy the line
	echo $line
	else
		# get gene name
		input_gene=$(echo $line | awk '{print $1}')
		# get exons of input gene 				  | convert to bed 								  | sort			   | remove duplicate lines			> store as temporary file for this gene
		grep $input_gene $REF$REFbase"_EXONS.gtf" | awk '{print $1"\t"$4-1"\t"$5"\texon\t.\t"$7}' | bedtools sort -i - | awk '!seen[$2$3]++ {print $0}' > "tmp."$SNa"_"$input_gene".all_exons.bed"
		# create bedfile with chromosomal location of exon-intron junction
		echo $line | awk '{print $2"\t"$11-1"\t"$11"\t"$1"\t.\t"$3}' > "tmp."$SNa"_out1_breakpoint.bed"
		# create list with genomic distances of exons to chromosomal breakpoint -strand specific -print distance to (a) -show top 100
		bedtools closest -s -D a -k 100 -a "tmp."$SNa"_out1_breakpoint.bed" -b "tmp."$SNa"_"$input_gene".all_exons.bed" > "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv"
		# determine genomic region of exon that constitutes the chromosomal breakpoint (this assumes that the correct exon comes out on top after running bedtools closest) | convert to "region" syntax in samtools
		region1=$(head -n1 "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{print $7":"$8"-"$9}')
		# check if $region1 is empty
		if [[ -z $region1 ]]
		then
			echo $line
		else
			# make sure collecting file does not yet exist
			rm -f "tmp."$SNa"_out5_TEreads"
			rm -f "tmp."$SNa"_out6_genereads"
			for ((snum=1; snum <= NumSam ; snum++))
			do
				# make sure collecting file does not yet exist
				rm -f "tmp."$SNa"_S"$snum"_out3_TEreads"
				rm -f "tmp."$SNa"_S"$snum"_out4_genereads"
				for ((l=1; l <= NumLan ; l++))
				do
					# check if STAR output bam has already been indexed, index if not
					if [[ -f $path_to_TEchimPART2_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam.bai" ]] ; then : ; else samtools index $path_to_TEchimPART2_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam" ; fi
					# get number of reads inside region 1 where mate maps onto transposon
					samtools view $path_to_TEchimPART2_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v te="$(echo $line | awk '{print $5}')" '{ if($7 ~ te) {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_S"$snum"_out3_TEreads"
					# for the next step, both the strandedness of the gene and the type of fragment (GENE-TE or TE-GENE) determines the necessary steps
					if [[ $(echo $line | awk '{print $3}') == "+" ]]
					then							
						if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]
						then
							# determine the nearest intron-exon junction beyond the breakpoint. if the gene is on the positive strand and the fragment is GENE-TE, then the next non-TE exon is located downstream (the smallest positive value in the bedtools closest output)
							next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13>0) {print$8}}' | head -n1)
							# first, get all reads that map in region1																	| next, extract (and count) reads where mate maps beyond the next exon junction.
							samtools view $path_to_TEchimPART2_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && $8>=ne && $6 == "60M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_S"$snum"_out4_genereads"
						elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]
						then
							# if fragment is TE-GENE, then the "next" intron-exon junction is upstream of the breakpoint ( the smallest NEGATIVE value in the bedtools closest output)
							next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13<0) {print$9}}' | head -n1)
							# first, get all reads that map in region1																	| next, extract (and count) reads where mate maps before the next exon junction.
							samtools view $path_to_TEchimPART2_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && $8+60<=ne && $6 == "60M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_S"$snum"_out4_genereads"
						fi
					elif [[ $(echo $line | awk '{print $3}') == "-" ]]
					then 
						if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]
						then
							# if the gene is on the negative strand	and the fragment is GENE-TE, then the next non-TE exon is upstream of breakpoint (BUT: TAKE smallest POSITIVE value in bedtools closest output)
							next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13>0) {print$9}}' | head -n1)
							# first, get all reads that map in region1																	| next, extract (and count) reads where mate maps upstream to the next exon junction.
							samtools view $path_to_TEchimPART2_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && ne!= "" && $8+60<=ne && $6 == "60M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_S"$snum"_out4_genereads"
						elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]
						then
							next_exon=$(cat "tmp."$SNa"_"$input_gene"_out2_exon_breakpoint_distances.tsv" | awk '{if ($13<0) {print$8}}' | head -n1)
							samtools view $path_to_TEchimPART2_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && ne!= "" && $8>=ne && $6 == "60M") {print $0}}' | wc -l | awk '{print $1}' >> "tmp."$SNa"_S"$snum"_out4_genereads"	
						fi
					fi
				done
				awk '{s+=$1} END {print s}' "tmp."$SNa"_S"$snum"_out3_TEreads" >> "tmp."$SNa"_out5_TEreads"
				awk '{s+=$1} END {print s}' "tmp."$SNa"_S"$snum"_out4_genereads" >> "tmp."$SNa"_out6_genereads"
			done
			col1=$(paste -sd "|" "tmp."$SNa"_out5_TEreads")
			col2=$(paste -sd "|" "tmp."$SNa"_out6_genereads")
			echo $line | awk -v c1="$col1" -v c2="$col2" '{print $0"\tTE:"c1"\tGene:"c2}'
		fi
		rm -f "tmp."$SNa*
	fi
done < $PART2_output > $SNa"_out20_withGENEreads.tsv"
