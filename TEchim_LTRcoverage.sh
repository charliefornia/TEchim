#!/bin/bash

################################################################################
# TITLE: TEchim - detection of TE-gene chimera in gDNA data
# VERSION: 0.1.2 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 42/05/2019 (dd/mm/yyyy)
# DESCRIPTION:
################################################################################

################################################################################
# REQUIREMENTS:
# - samtools
################################################################################

################################################################################
################################################################################
# set parameters
wd=/Users/koalcan/Dropbox/CloudDesktop/TEchim_cloud/DEV/
SNa=2018MARCH_TEchim
REFpath=/Users/koalcan/Documents/REF_2019FEB/for_TEchim/
REFbase=dmel625
TElist=TEs.fa
path_to_PART1_output=/Users/koalcan/Documents/2019MAY_TEscoex_longRNA/fromHPC/
NumSam=6
NumLan=2
################################################################################
################################################################################

# change to wd
cd $wd

if [[ -f $REFpath$TElist".list_of_LTR_TE_names.txt" ]]; then :; else awk '{if ($2 == "LTR" && $1 !~ "_LTR") {gsub(/>/,""); print $1}}' $REFpath$TElist > $REFpath$TElist".list_of_LTR_TE_names.txt"; fi

for ((snum=1; snum <= NumSam ; snum++))
do
	for ((l=1; l <= NumLan ; l++))
	do
	if [[ -f $SNa".LTRinput.S"$snum"_L"$l".sam" ]]; then :; else samtools view $path_to_PART1_output$SNa"_S"$snum"_L"$l"_out4_Aligned.sortedByCoord.out.bam" | awk '{if ($3 ~ "TEchr_" && $7 != "=") print $0}' > $SNa".LTRinput.S"$snum"_L"$l".sam"; fi
	done		
done

while read line
do
	rm -f "tmp.output.collect_"$line
	for ((snum=1; snum <= NumSam ; snum++))
	do	
		rm -f "tmp.output.S"$snum
		for ((l=1; l <= NumLan ; l++))
		do
		awk -v te="$line" '{if ($3 == "TEchr_"te"_LTR" && $7 == "TEchr_"te) counta++; else if ($3 == "TEchr_"te"_LTR" && $7 != "TEchr_"te && $7 != "=") countb++} END {print countb"\t"counta}' $SNa".LTRinput.S"$snum"_L"$l".sam" >> "tmp.output.S"$snum
		done
		awk '{suma+=$1; sumb+=$2;} END {print suma"\t"sumb;}' "tmp.output.S"$snum >> "tmp.output.collect_"$line
		rm -f "tmp.output.S"$snum		
	done
	a=$(cat "tmp.output.collect_"$line | awk '{print $1}' | paste -sd "|" -)
	b=$(cat "tmp.output.collect_"$line | awk '{print $2}' | paste -sd "|" -)
	echo $line $a $b
	rm -f "tmp.output.collect_"$line
done < $REFpath$TElist".list_of_LTR_TE_names.txt" > $SNa".Samples1to"$NumSam".LTR_TEs.proportion_of_TEvsGENE.tsv"
