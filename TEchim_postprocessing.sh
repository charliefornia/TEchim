#!/bin/bash

################################################################################
# TITLE: TEchim - IGE
# VERSION: 0.2.2 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 19/06/2019 (dd/mm/yyyy)
# DESCRIPTION: Postprocessing 
################################################################################

################################################################################
################################################################################
# set parameters
wd=$(pwd)						# working directory (no trailing "/")
path_to_PART1_output=$wd"/"		# with trailing "/"
nc=1							# number of cores
paralogs=/PATH/TO/paralogs.gff 
################################################################################
################################################################################
# functions:

get_variables()
{
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
}

#### POST PROCESSING of TEchim


# 1.) Filter reads


################################################################################
################################################################################
# execute:

# change to wd
cd $wd
get_variables
cd $wd"/IGE_COLLECTION_"$SNa"/"
input_letters=$(printf "a\tb\tc\td\te\tf\tg\th\ti\tj")	
for IGEgroup in $input_letters
do
	(while read line; do a=$(echo $line | awk '{print$1}'); b=$(echo $line | awk '{gsub ("@","\t"); print $5}'); c=$(grep $a $paralogs | grep $b | wc -l | awk '{print $1}'); if [[ $c = "0" ]]; then echo $line | awk '{gsub(" ","\t"); print $0}'; fi; done < $IGEgroup"_"$SNa"_IGEref_chimericreads_final.FILTERED.tsv" > $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP2"
	# remove cases where the hit is within +/- 1000nt of IGE
	awk '{ orig=$0; gsub("@","\t",$5); a=$5; gsub("-","\t",a); $0=a; gsub("\\(","\t",$4); gsub(":","\t"); chr=$3; start=$4-1000; end=$5+1000; $0=orig; if($2!=chr) print $0; else if (start > $4 && end > $4) print $0}' $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP2" > $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP3"	
	awk '{print $2"\t"$4-21"\t"$4+20"\tLINE"NR"\t.\t"$3}' $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP3" > $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP4.bed"
	list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
	for SNo in $list_of_SNo
	do
		list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
		for LNo in $list_of_LNo
		do
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" -bed $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP4.bed" -s  | awk '{print $7}' > "tmp."$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_insertionsites.coverage_perSample.tsv" &
		done
	done
	wait
	paste "tmp."$IGEgroup"_"$SNa"_"*"_insertionsites.coverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$IGEgroup"_"$SNa"_IGEcoverage_averages.tsv"
	# combine averages with chimeric output
	paste $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP3" "tmp."$IGEgroup"_"$SNa"_IGEcoverage_averages.tsv" > $IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP4.tsv") &
done
wait


cd $wd
awk '{print $2"\t"$4-21"\t"$4+20"\tLINE"NR"\t.\t"$3}' $SNa"_TE_chimericreads_final.withGENES.ABOVE1.tsv" > $SNa"_TE_chimericreads_final.STEP4.bed"
list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
for SNo in $list_of_SNo
do
	list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
	for LNo in $list_of_LNo
	do
		bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" -bed $SNa"_TE_chimericreads_final.STEP4.bed" -s  | awk '{print $7}' > "tmp."$SNa"_"$SNo"_"$LNo"_TE_insertionsites.coverage_perSample.tsv" &
	done
done
wait
paste "tmp."$SNa"_"*"_TE_insertionsites.coverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$SNa"_TEcoverage_averages.tsv"
# combine averages with chimeric output
paste $SNa"_TE_chimericreads_final.withGENES.ABOVE1.tsv" "tmp."$SNa"_TEcoverage_averages.tsv" > $SNa"_TEref_chimericreads_final.STEP4.tsv"




