#!/bin/bash

################################################################################
# TITLE: TEchim - quantify chimeras using LTR-spanning reads
# VERSION: 0.2.0 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 13/06/2019 (dd/mm/yyyy)
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
REFpath=						# if left empty then REFpath from PART1 is used
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

# assign TElist parameter
if [ -e $REFpath"REFERENCE_TElist" ]; then
	TElist=$(cat $REFpath"REFERENCE_TElist")
else
	echo " #### ERROR: reference path is corrupt - no file named REFERENCE_TElist"
	exit
fi

# create list of all LTR TEs
if [[ -f $REFpath$TElist".list_of_LTR_TE_names.txt" ]]; then :; else awk '{if ($2 == "LTR" && $1 !~ "_LTR") {gsub(/>/,""); print $1}}' $REFpath$TElist > $REFpath$TElist".list_of_LTR_TE_names.txt"; fi

####################################################

list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
for SNo in $list_of_snum
do
	list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
	for LNo in $list_of_lanes
	do
		if [[ -f "tmp."$SNa".LTRinput."$SNo"_"$LNo".sam" ]]; then :; else samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" | awk '{if ($3 ~ "TEchr_" && $7 != "=") print $0}' > "tmp."$SNa".LTRinput."$SNo"_"$LNo".sam"; fi
	done		
done

# write header
echo TE$'\t'LTR-gene$'\t'LTR-TE > $SNa".LTR_TEs.proportion_of_TEvsGENE.tsv"

while read line
do
	rm -f "tmp.output.collect_"$line
	list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
	for SNo in $list_of_snum
	do
		rm -f "tmp.output."$SNo
		list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
		for LNo in $list_of_lanes
		do
			awk -v te="$line" '{if ($3 == "TEchr_"te"_LTR" && $7 == "TEchr_"te) counta++; else if ($3 == "TEchr_"te"_LTR" && $7 != "TEchr_"te && $7 != "=") countb++} END {print countb"\t"counta}' "tmp."$SNa".LTRinput."$SNo"_"$LNo".sam" >> "tmp.output."$SNo
		done
		awk '{suma+=$1; sumb+=$2;} END {print suma"\t"sumb;}' "tmp.output."$SNo >> "tmp.output.collect_"$line
		rm -f "tmp.output."$SNo		
	done
	a=$(cat "tmp.output.collect_"$line | awk '{print $1}' | paste -sd "|" -)
	b=$(cat "tmp.output.collect_"$line | awk '{print $2}' | paste -sd "|" -)
	echo $line $a $b
	rm -f "tmp.output.collect_"$line
done < $REFpath$TElist".list_of_LTR_TE_names.txt" >> $SNa".LTR_TEs.proportion_of_TEvsGENE.tsv"
rm "tmp."$SNa".LTRinput."*
