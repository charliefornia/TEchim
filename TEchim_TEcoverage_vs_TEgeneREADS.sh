#!/bin/bash

################################################################################
# TITLE: TEchim - quantify chimeras using LTR-spanning reads
# VERSION: 0.3.0 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 09/04/2020 (dd/mm/yyyy)
# DESCRIPTION: This tool measures the difference between LTR-TE vs. LTR-gene
# breakpoint-spanning contigs
################################################################################

################################################################################
# REQUIREMENTS:
# - samtools
################################################################################

################################################################################
################################################################################
# set parameters
wd=$(pwd)
path_to_PART1_output=$wd"/"
################################################################################
################################################################################

# change to wd
cd $wd

# assign sample name
if [ -e $path_to_PART1_output"."*"_samplename" ]; then
	if [ $(cat $path_to_PART1_output"."*"_samplename" | wc -l | awk '{print $1}') = 1 ]; then
		SNa=$(cat $path_to_PART1_output"."*"_samplename")
	else
		echo " #### ERROR: path to output from PART1 is corrupt - more than one file named *_samplename"
	fi
else
	echo " #### ERROR: path to output from PART1 is corrupt - no file named *_samplename"
	exit
fi


cd $wd
list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
for SNo in $list_of_SNo
	do
	list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
	for LNo in $list_of_LNo
		do
		samtools view -f 3 -h $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" | grep TEchr | samtools view -Sb | bedtools genomecov -ibam - -d | cat - <(echo "end") | awk '$1==last {sum+=$3; count+=1; next} NR>1 {print " ",sum/count;} {last=$1; printf "%s",$1,sum; sum=0; count=0}' | head -n-1 > "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEaverage"
		samtools view -f 1 -F 2 -h $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" | awk '{if ($3 ~ "TEchr" && $7 !~ "TEchr" && $7 != "=") print $0; else if ($3 !~ "TEchr" && $7 ~ "TEchr") print $0 }' > "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEtogeneREADS"
		if [[ -f "tmp.TEcoverage."$SNa"_TEnames" ]] ; then : ; else awk '{print $1}' "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEaverage" > "tmp.TEcoverage."$SNa"_TEnames" ; fi
		while read line
		do
			awk -v te="$(echo $line | awk '{print $1}')" '{if ($3 == te) print $1; else if ($7 == te) print $1}' "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEtogeneREADS" | wc -l | paste <(echo $line | awk '{print $2}') -
		done < "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEaverage" > "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEaveragesANDgeneREADS"
		rm "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEaverage"
		rm "tmp.TEcoverage."$SNa"_"$SNo"_"$LNo"_TEtogeneREADS"
		done
	done
paste "tmp.TEcoverage."$SNa"_TEnames" "tmp.TEcoverage."$SNa"_"*"_TEaveragesANDgeneREADS" > $SNa"_TEcoverage_and_TEgeneREADS.tsv"
	
