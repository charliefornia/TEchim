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

# create list of all LTR TEs (the sed command removes "_I" from some TE names - this is so the TE-base name matches the name of its LTR)
awk '{if ($2 == "LTR" && $1 !~ "_LTR") {gsub(/>/,""); print $1}}' $REFpath$TElist | sed -e 's/\(_I\)*$//g' > $TElist".list_of_LTR_TE_names.txt"

####################################################

list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
for SNo in $list_of_SNo
do
	list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
	for LNo in $list_of_LNo
	do
		# extract BAM header
		samtools view -H $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" > "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam"
		samtools view -H $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" > "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam"
		
		# extract reads where FIRST read maps to LTR region
		samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" | awk '{if ($3 ~ "TEchr_" && $3 ~ "_LTR" && $7 != "=") print $0}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam"
		# extract reads where SECOND read maps to LTR region
		samtools view $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" | awk '{if ($7 ~ "TEchr_" && $7 ~ "_LTR") print $0}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam"
			
		# Select pairs where both reads map to positive strand
		samtools view -f 65 -F 48 "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam" | awk '{if ($7 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_B.sam"
		samtools view -f 65 -F 48 "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam" | awk '{if ($3 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_A.sam"
		
		# Select pairs where read1 maps to positive, and read2 to negative
		samtools view -f 97 -F 16 "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam" | awk '{if ($7 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_B.sam"
		samtools view -f 97 -F 16 "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam" | awk '{if ($3 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_A.sam"
		samtools view -f 97 -F 16 "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam" | awk '{if ($3 == $7"_LTR" ) {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_D.sam"
		samtools view -f 97 -F 16 "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam" | awk '{if ($7 == $3"_LTR" ) {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_C.sam"
		
		# Select pairs where read1 maps to negative, and read2 to positive
		samtools view -f 81 -F 32 "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam" | awk '{if ($7 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_A.sam"
		samtools view -f 81 -F 32 "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam" | awk '{if ($3 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_B.sam"
		samtools view -f 81 -F 32 "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam" | awk '{if ($3 == $7"_LTR" ) {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_C.sam"
		samtools view -f 81 -F 32 "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam" | awk '{if ($7 == $3"_LTR" ) {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_D.sam"		
		
		# Select pairs where both reads map to negative strand
		samtools view -f 113 "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam" | awk '{if ($7 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_A.sam"
		samtools view -f 113 "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam" | awk '{if ($3 !~ "TEchr_") {print $0}}' >> "tmp."$SNa".LTRinput."$SNo"_"$LNo".COLLECT_B.sam"
		
		rm "tmp."$SNa".LTRinput."$SNo"_"$LNo".FIRSTREADISLTR.sam"
		rm "tmp."$SNa".LTRinput."$SNo"_"$LNo".SECONDREADISLTR.sam"
	
	done		
done

wait

# write header
echo TE$'\t'GENE-LTR$'\t'LTR-GENE$'\t'TE-LTR$'\t'LTR-TE > $SNa".LTR_TEs.proportion_of_TEvsGENE.tsv"

rm -f "tmp."$SNa".output.collect."*".txt"
	
while read line
do		
	list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
	for SNo in $list_of_SNo
	do
		cat "tmp."$SNa".LTRinput."$SNo"_"*".COLLECT_A.sam" > "tmp."$SNa".LTRTEs."$SNo".COLLECT_A.sam"
		cat "tmp."$SNa".LTRinput."$SNo"_"*".COLLECT_B.sam" > "tmp."$SNa".LTRTEs."$SNo".COLLECT_B.sam"
		cat "tmp."$SNa".LTRinput."$SNo"_"*".COLLECT_C.sam" > "tmp."$SNa".LTRTEs."$SNo".COLLECT_C.sam"
		cat "tmp."$SNa".LTRinput."$SNo"_"*".COLLECT_D.sam" > "tmp."$SNa".LTRTEs."$SNo".COLLECT_D.sam"
		
		awk -v te=$line '{if ($3 ~ te || $7 ~ te) count++1} END {printf count"|"}' "tmp."$SNa".LTRTEs."$SNo".COLLECT_A.sam" >> "tmp."$SNa".output.collect.A.txt"
		awk -v te=$line '{if ($3 ~ te || $7 ~ te) count++1} END {printf count"|"}' "tmp."$SNa".LTRTEs."$SNo".COLLECT_B.sam" >> "tmp."$SNa".output.collect.B.txt"
		awk -v te=$line '{if ($3 ~ te || $7 ~ te) count++1} END {printf count"|"}' "tmp."$SNa".LTRTEs."$SNo".COLLECT_C.sam" >> "tmp."$SNa".output.collect.C.txt"
		awk -v te=$line '{if ($3 ~ te || $7 ~ te) count++1} END {printf count"|"}' "tmp."$SNa".LTRTEs."$SNo".COLLECT_D.sam" >> "tmp."$SNa".output.collect.D.txt"
	done
	echo $line$'\t'$(head -n1 "tmp."$SNa".output.collect.A.txt")$'\t'$(head -n1 "tmp."$SNa".output.collect.B.txt")$'\t'$(head -n1 "tmp."$SNa".output.collect.C.txt")$'\t'$(head -n1 "tmp."$SNa".output.collect.D.txt")
	rm -f "tmp."$SNa".LTRTEs."*
	rm -f "tmp."$SNa".output.collect."*".txt"	
done < $TElist".list_of_LTR_TE_names.txt" >> $SNa".LTR_TEs.proportion_of_TEvsGENE.tsv"

rm "tmp."$SNa".LTRinput."*
