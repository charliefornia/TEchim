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

# for IGEs
cd $wd
get_variables
cd $wd"/IGE_COLLECTION_"$SNa"/"
input_letters=$(printf "a\tb\tc\td\te\tf\tg\th\ti\tj")	
for IGEgroup in $input_letters
do
	# remove mitochondrial hits and hits in less than 2 replicates
	(grep -v mito $IGEgroup"_"$SNa"_IGEref_chimericreads_final.FILTERED.tsv" | awk '{ if ($8 > 1) {print $0}}' > "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP1"
	while read line; do a=$(echo $line | awk '{print$1}'); b=$(echo $line | awk '{gsub ("@","\t"); print $5}'); c=$(grep $a $paralogs | grep $b | wc -l | awk '{print $1}'); if [[ $c = "0" ]]; then echo $line | awk '{gsub(" ","\t"); print $0}'; fi; done < "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP1" > "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP2"
	# remove cases where the hit is within +/- 5kb of IGE
	awk '{ orig=$0; gsub("@","\t",$5); a=$5; gsub("-","\t",a); $0=a; gsub("\\(","\t",$4); gsub(":","\t"); chr=$2; start=$3-5000; end=$4+5000; $0=orig; if($2!=chr) print $0; else if (start > $4 || end < $4) print $0}' "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP2" > "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP3"	
	awk '{print $2"\t"$4-21"\t"$4+20"\tLINE"NR"\t.\t"$3}' "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP3" > "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP4.bed"
	list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
	for SNo in $list_of_SNo
	do
		list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
		for LNo in $list_of_LNo
		do
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" -bed "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP4.bed" -s  | awk '{print $7}' > "pp.tmp."$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_insertionsites.coverage_perSample.tsv" &
		done
	done
	wait
	rm "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP4.bed"
	paste "pp.tmp."$IGEgroup"_"$SNa"_"*"_insertionsites.coverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "pp.tmp."$IGEgroup"_"$SNa"_IGEcoverage_averages.tsv"
	rm "pp.tmp."$IGEgroup"_"$SNa"_"*"_insertionsites.coverage_perSample.tsv"
	# combine averages with chimeric output
	paste "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP3" "pp.tmp."$IGEgroup"_"$SNa"_IGEcoverage_averages.tsv" > $IGEgroup"_"$SNa"_IGEref_chimericreads_final.FILTERED.ABOVE1.withAverages.tsv"
	rm "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP1"
	rm "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP2"
	rm "pp.tmp."$IGEgroup"_"$SNa"_IGEref_chimericreads_final.STEP3"
	rm "pp.tmp."$IGEgroup"_"$SNa"_IGEcoverage_averages.tsv"
	) &
done
wait


# for TEs
cd $wd
awk '{if ($8 > 1) {print "LINE"NR"\t"$0}}' $SNa"_TE_chimericreads_final_withGENEreads.tsv" > "pp.tmp."$SNa"_TE_chimericreads_final.ABOVE1.STEP3.tsv"
awk '{print $3"\t"$5-21"\t"$5+20"\t"$1"\t.\t"$4}' "pp.tmp."$SNa"_TE_chimericreads_final.ABOVE1.STEP3.tsv" > "pp.tmp."$SNa"_TE_chimericreads_final.STEP4.bed"
list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
for SNo in $list_of_SNo
do
	list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
	for LNo in $list_of_LNo
	do
		bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" -bed "pp.tmp."$SNa"_TE_chimericreads_final.STEP4.bed" -s  | awk '{print $7}' > "pp.tmp."$SNa"_"$SNo"_"$LNo"_TE_insertionsites.coverage_perSample.tsv" &
	done
done
wait
rm "pp.tmp."$SNa"_TE_chimericreads_final.STEP4.bed"
paste "pp.tmp."$SNa"_"*"_TE_insertionsites.coverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "pp.tmp."$SNa"_TEcoverage_averages.tsv"
rm "pp.tmp."$SNa"_"*"_TE_insertionsites.coverage_perSample.tsv"
paste "pp.tmp."$SNa"_TE_chimericreads_final.ABOVE1.STEP3.tsv" "pp.tmp."$SNa"_TEcoverage_averages.tsv" | cut -f 2- > $SNa"_TE_chimericreads_final.withGENES.ABOVE1.withAverages.tsv"
rm "pp.tmp."$SNa"_TEcoverage_averages.tsv"
rm "pp.tmp."$SNa"_TE_chimericreads_final.ABOVE1.STEP3.tsv"
cat <(cat $SNa"_TE_chimericreads_final.withGENES.ABOVE1.withAverages.tsv") <(awk '{if ($8 <= 1) {print $0}}' $SNa"_TE_chimericreads_final_withGENEreads.tsv") > $SNa"_TE_chimericreads_final.withGENES.ALL.withAverages.tsv"








