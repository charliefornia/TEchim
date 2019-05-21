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
SNa=2018MARCH_TEchim
SNo=1
LNo=1
REF=~/Dropbox/CloudDesktop/REF_cloud/
REFbase=dmel625_v04
################################################################################

# check if STAR output bam has already been indexed, index if not
if [[ -f $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_out4_Aligned.sortedByCoord.out.bam.bai" ]] ; then : ; else samtools index $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_out4_Aligned.sortedByCoord.out.bam" ; fi
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
		input_gene=$(echo $line | awk '{print $1}')
		grep $input_gene $REF$REFbase"_EXONS.gtf" | awk '{print $1"\t"$4-1"\t"$5"\texon\t.\t"$7}' | bedtools sort -i - | awk '!seen[$2$3]++ {print $0}' > "tmp."$input_gene".all_exons.bed"
		echo $line | awk '{print $2"\t"$11-1"\t"$11"\t"$1"\t.\t"$3}' > "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out1_breakpoint.bed"
		bedtools closest -s -D a -a "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out1_breakpoint.bed" -b "tmp."$input_gene".all_exons.bed" -k 100 > "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out2_exons_near_breakpoint.bed"
		region1=$(head -n1 "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out2_exons_near_breakpoint.bed" | awk '{print $7":"$8"-"$9}')
		out_te=$(samtools view $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v te="$(echo $line | awk '{print $5}')" '{ if($7 ~ te) {print $0}}' | wc -l | awk '{print $1}')
		if [[ $(echo $line | awk '{print $3}') == "+" ]]
		then							
			if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]
			then
				next_exon=$(cat "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out2_exons_near_breakpoint.bed" | awk '{if ($13>0) {print$8}}' | head -n1)
				samtools view $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && $8>=ne && $6 == "60M") {print $0}}' | wc -l | awk -v line="$line" -v sno="$SNo" -v lno="$LNo" -v nte="$out_te" '{print line"\tS"sno"_L"lno":"nte"/"$1}'
			elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]
			then
				next_exon=$(cat "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out2_exons_near_breakpoint.bed" | awk '{if ($13<0) {print$9}}' | head -n1)
				samtools view $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && $8<=ne && $6 == "60M") {print $0}}' | wc -l | awk -v line="$line" -v sno="$SNo" -v lno="$LNo" -v nte="$out_te" '{print line"\tS"sno"_L"lno":"nte"/"$1}'
			fi
		elif [[ $(echo $line | awk '{print $3}') == "-" ]]
		then 
			if [[ $(echo $line | awk '{print $6}') == "GENE-TE" ]]
			then
				next_exon=$(cat "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out2_exons_near_breakpoint.bed" | awk '{if ($13>0) {print$9}}' | head -n1)
				samtools view $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && ne!= "" && $8<=ne && $6 == "60M") {print $0}}' | wc -l | awk -v line="$line" -v sno="$SNo" -v lno="$LNo" -v nte="$out_te" '{print line"\tS"sno"_L"lno":"nte"/"$1}'
			elif [[ $(echo $line | awk '{print $6}') == "TE-GENE" ]]
			then
				next_exon=$(cat "tmp."$input_gene"."$SNa"_S"$SNo"_L"$LNo"_out2_exons_near_breakpoint.bed" | awk '{if ($13<0) {print$8}}' | head -n1)
				samtools view $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_out4_Aligned.sortedByCoord.out.bam" $region1 | awk -v ne="$next_exon" '{ if($7 == "=" && ne!= "" && $8>=ne && $6 == "60M") {print $0}}' | wc -l | awk -v line="$line" -v sno="$SNo" -v lno="$LNo" -v nte="$out_te" '{print line"\tS"sno"_L"lno":"nte"/"$1}'	
			fi
		fi
		rm -f "tmp."$input_gene*
	fi
done < $path_to_TEchimPART2_output$SNa"_S"$SNo"_L"$LNo"_chimericreads_final.tsv" > $SNa"_S"$SNo"_L"$LNo"_out20_withGENEreads.tsv"



