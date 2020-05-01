#!/bin/bash

################################################################################
# TITLE: TEchim - PART 2
# VERSION: 0.4.0 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 18/03/2020 (dd/mm/yyyy)
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
	
	if [ -e $path_to_PART1_output"."$SNa"_maxfraglength" ]; then
		MaxFragLength=$(sort -nr $path_to_PART1_output"."$SNa"_maxfraglength" | head -n1)
	else
		echo " #### ERROR: path to output from PART1 is corrupt - no file named .""$SNa""_maxfraglength"
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
	echo "=====================" > $wd"/"$SNa"_PART2_"$logname".log"
	echo "|| TEchim - PART 2 || " >> $wd"/"$SNa"_PART2_"$logname".log"
	echo "=====================" >> $wd"/"$SNa"_PART2_"$logname".log"
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
	cat $path_to_PART1_output$SNa*"/"$SNa*"_out13_breakpoints.bed" | bedtools sort -i - | tr '|' '\t' > $SNa"_out02_sepparated.tsv" 
	# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
	# and create column with BASE readname (no :A or :B)
	awk 'BEGIN {OFS = "\t"} {a = $5 ; gsub(/_LTR/,"",a) ; c = substr($4, 1, length($4)-2) ; print $0"\t"a"\t"c}' < $SNa"_out02_sepparated.tsv" > $SNa"_out03pre_TEbase.tsv" && rm $SNa"_out02_sepparated.tsv"
	wc_predup=$(wc -l $SNa"_out03pre_TEbase.tsv" | awk '{print $1}')
	# remove duplicate reads (where both :A and :B version of the same reads were picked up
	awk '!seen[$15]++' $SNa"_out03pre_TEbase.tsv" > $SNa"_out03_TEbase.tsv" && rm $SNa"_out03pre_TEbase.tsv"
	rm $SNa"_in10_combined.sorted.bed"
	echo " ---- $(wc -l $SNa"_out03_TEbase.tsv" | awk '{print $1}') TE reads ($wc_predup before removing duplicates)" >> $wd"/"$SNa"_PART2_"$logname".log"
	echo " <-- done processing PART1 output for TE at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
}

combine_hits_of_each_TE()
{
	cd $wd
	echo " --> start combining TE reads at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
	# create list of all TEs in dataset
	cut -f14 $1 | sort | uniq > $SNa"_out03a_uniqueTEs.tsv"
	# to create collection file, first make sure this file does not yet exist
	rm -f $SNa"_out12_OUTPUT.tsv"
	# loop through all TEs in the dataset
	while read TE
	do
		# grep TEs from combined data set
		awk -v t="$TE" '{if ($14 == t) print $0;}' $1 > "tmp."$SNa"_"$TE"_out04.tsv"
		if [ -s "tmp."$SNa"_"$TE"_out04.tsv" ]; then	
			###### SPLIT into upstream and downstream
			awk '{if ($7 == "upstream") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_upstream.tsv"
			# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
			if [ -s "tmp."$SNa"_"$TE"_out05_upstream.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_upstream.tsv" -c 5,12,13,8,15,3,6,7,8,10,5 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,collapse,collapse,distinct -d 20 > "tmp."$SNa"_"$TE"_out06_upstream.tsv"
			fi
			awk '{if ($7 == "downstream") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_downstream.tsv"
			if [ -s "tmp."$SNa"_"$TE"_out05_downstream.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_downstream.tsv" -c 5,12,13,8,15,3,6,7,8,10,5 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,collapse,collapse,distinct -d 20 > "tmp."$SNa"_"$TE"_out06_downstream.tsv"
			fi
			cat "tmp."$SNa"_"$TE"_out06_"*".tsv" | bedtools sort -i - | awk 'BEGIN {OFS = "\t"} {print "gDNA",$1,".",$9,$4,$11,$10,$7,$8,"gDNA",$12,$13,$14 }' >> $SNa"_out12_OUTPUT.tsv"
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

add_features()
{
	cd $wd
	echo " --> start checking for Features at breakpoint ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
	# now, every line is one incident of a downstream chimera
	while read line
	do
		# convert line to bed format
		echo $line | awk '{ print $2"\t"$4-1"\t"$4"\t"$5"|"$1"|"$6"\t.\t"$3}' > "tmp."$SNa"_out13.bed"
		d=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa | awk '{print $10"("$12")"}' | paste -sd ";" -)
		echo $line | awk -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"d}'
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
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$13,$14}' $1 > $SNa"_newcola"
	# append pooled TE breakpoints to final output
	paste $SNa"_newcola" $SNa"_newcolb" $SNa"_newcolc" > $SNa"_TE_chimericreads_final.tsv"
	rm $SNa"_newcola" $SNa"_newcolb" $SNa"_newcolc"
	rm $1
	echo " <-- done tyding up TE breakpoints at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
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
add_features $wd"/"$SNa"_out12a_OUTPUT.tsv"
split_TE_breakpoints $wd"/"$SNa"_out15.tsv"



echo " <-- all done at ... $(date)" >> $wd"/"$SNa"_PART2_"$logname".log"
echo "================================" >> $wd"/"$SNa"_PART2_"$logname".log"
