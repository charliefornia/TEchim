#!/bin/bash

################################################################################
# TITLE: TEchim - IGE
# VERSION: 0.4.0 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 18/03/2020 (dd/mm/yyyy)
# DESCRIPTION: This script generates 10 IGE subsamples 
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
nc=1							# number of cores
REFpath=						# if left empty then REFpath from PART1 is used
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
	echo "==================" > $wd"/"$SNa"_IGE_"$logname".log"
	echo "|| TEchim - IGE || " >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "==================" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "Parameters:" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "Working directory:" "$wd" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "Location of PART1 output:" "$path_to_PART1_output" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "Sample name:" "$SNa" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "Reference files:" "$REFpath$REFbase" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "--------------------------------" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo " --> starting at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo "--------------------------------" >> $wd"/"$SNa"_IGE_"$logname".log"
}

split_CDS()
{
	cd $wd
	# if clause to only create this directory once
	if [ ! -d "IGE_COLLECTION_"$SNa ]; then
		mkdir "IGE_COLLECTION_"$SNa
		cd "IGE_COLLECTION_"$SNa
		# randomly subsample filtered CDS (here: take 10%)
		tenpc_CdS=$(wc -l $REFpath$REFbase".CDS_for_IGE.bed" | awk '{OFMT = "%.0f"; print $1/10}')
		# split into ten random groups
		sort -R $REFpath$REFbase".CDS_for_IGE.bed" | split -l $tenpc_CdS -a 1 - $SNa"_inputGENEs_"	
	fi
}

calculate_TE_coverage()
{
	cd $wd"/IGE_COLLECTION_"$SNa
	if [[ ! -f $REFpath"TEs.POSITIVE.fa.bed" ]]; then awk '{if ($6 == "+") print $0}' $REFpath"TEs.fa.bed" > $REFpath"TEs.POSITIVE.fa.bed"; fi
	if [[ ! -f $REFpath"TEs.NEGATIVE.fa.bed" ]]; then awk '{if ($6 == "-") print $0}' $REFpath"TEs.fa.bed" > $REFpath"TEs.NEGATIVE.fa.bed"; fi
	# create list with all sample numbers in PART1 output
	list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
	for SNo in $list_of_SNo
	do
		# create list with lane numbers for each sample in PART1 output
		list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
		for LNo in $list_of_LNo
		do
			# check if .bam file is indexed, if not then index it
			if [[ ! -f $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam.bai" ]]; then samtools index $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam"; fi
			# create bam mapping to positive strand
			# POSITIVE_part1: -f 65 -> 1: read paired, 64: first in pair | -F 16 -> 16: read reverse strand
			# POSITIVE_part2: -f 145 -> 1: read paired, 16: read reverse strand, 128: second in pair
			samtools view -hb -f 65 -F 16 $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" > $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE_1.bam"
			samtools view -hb -f 145 $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" > $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE_2.bam"
			samtools merge -f $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE.bam" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE_1.bam" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE_2.bam"
			samtools index $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE.bam"
			# create bam mapping to negative strand
			# POSITIVE_part1: -f 129 -> 1: read paired, 128: second in pair | -F 16 -> 16: read reverse strand
			# POSITIVE_part2: -f 81 -> 1: read paired, 16: read reverse strand, 64: first in pair
			samtools view -hb -f 129 -F 16 $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" > $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE_1.bam"
			samtools view -hb -f 81 $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.bam" > $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE_2.bam"
			samtools merge -f $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE.bam" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE_1.bam" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE_2.bam"
			samtools index $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE.bam"
			rm $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/tmp."$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out."*
			# quantify coverage of TEs in each sample/lane - this line will print the number
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE.bam" -bed $REFpath"TEs.POSITIVE.fa.bed" | awk '{print $7}' > "tmp."$SNa"_"$SNo"_"$LNo"_out31_TEcoverage_perSample.tsv"
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE.bam" -bed $REFpath"TEs.NEGATIVE.fa.bed" | awk '{print $7}' >> "tmp."$SNa"_"$SNo"_"$LNo"_out31_TEcoverage_perSample.tsv"
		done
	done
	cat $REFpath"TEs.POSITIVE.fa.bed" $REFpath"TEs.NEGATIVE.fa.bed" | awk '{print $4}' > "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv"
	# combine all TE coverage numbers					  | calculate average, round to full integer
	paste "tmp."$SNa"_"*"_out31_TEcoverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$SNa"_out31b_TEcoverage_averages.tsv"
	# combine averages with TE names
	paste "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv" "tmp."$SNa"_out31b_TEcoverage_averages.tsv" | sort -k2n > $SNa"_TEcoverage_averages.tsv"
	rm "tmp."$SNa*
}

find_matching_IGEs()
{
	cd $wd"/IGE_COLLECTION_"$SNa
	echo " --> start finding matching IGE_$IGEgroup at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
	if [[ ! -f $SNa"_inputGENEs_POSITIVE_"$IGEgroup ]]; then awk '{if ($6 == "+") print $0}' $SNa"_inputGENEs_"$IGEgroup > $SNa"_inputGENEs_POSITIVE_"$IGEgroup; fi
	if [[ ! -f $SNa"_inputGENEs_NEGATIVE_"$IGEgroup ]]; then awk '{if ($6 == "-") print $0}' $SNa"_inputGENEs_"$IGEgroup > $SNa"_inputGENEs_NEGATIVE_"$IGEgroup; fi
	# assess how many Samples and Lanes exist - run bedtools mutlicov on each sample/lane
	list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
	for SNo in $list_of_SNo
	do
		list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
		for LNo in $list_of_LNo
		do
			# quantify coverage of CDS (subsample) in each sample/lane - this line will print the number
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.POSITIVE.bam" -bed $SNa"_inputGENEs_POSITIVE_"$IGEgroup | awk '{print $7}' > "tmp."$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_out31_GENEcoverage_perSample.tsv"
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo"_STAR/"$SNa"_"$SNo"_"$LNo"_out4_Aligned.sortedByCoord.out.NEGATIVE.bam" -bed $SNa"_inputGENEs_NEGATIVE_"$IGEgroup | awk '{print $7}' >> "tmp."$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_out31_GENEcoverage_perSample.tsv"
			done
	done
	cat $SNa"_inputGENEs_POSITIVE_"$IGEgroup $SNa"_inputGENEs_NEGATIVE_"$IGEgroup | awk '{print $4}' > "tmp."$IGEgroup"_"$SNa"_out31a_GENEcoverage_GENEnames.tsv"	
	# combine all GENE coverage numbers						| calculate average, round to full integer
	paste "tmp."$IGEgroup"_"$SNa"_"*"_out31_GENEcoverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$IGEgroup"_"$SNa"_out31b_GENEcoverage_averages.tsv"
	# combine averages with GENE names
	paste "tmp."$IGEgroup"_"$SNa"_out31a_GENEcoverage_GENEnames.tsv" "tmp."$IGEgroup"_"$SNa"_out31b_GENEcoverage_averages.tsv" | sort -k2n > $IGEgroup"_"$SNa"_GENEcoverage_averages.tsv"
	# loop through all TEs
	while read TEline
	do
		coverage=$(echo $TEline | awk '{print $2}')
		# find the gene that matches the expression level of each TE									  | print matching TE-GENE pairs
		awk -v v1="$coverage" '$2>v1 {print $1"\t"$2; exit}' $IGEgroup"_"$SNa"_GENEcoverage_averages.tsv" | awk -v te="$TEline" '{print te"\t"$0}'
	done < $wd"/IGE_COLLECTION_"$SNa"/"$SNa"_TEcoverage_averages.tsv" > "tmp."$IGEgroup"_"$SNa"_out35_TEmatchedCDS.tsv"
	# left-join TE-GENE pairs with original GENE bed (this is to get the strand of the according CDS's)																				 | print a bed file with each UNIQUE CDS (strand-specific!)
	join <(awk '{print $3"\t"$1"\t"$2"\t"$4}' "tmp."$IGEgroup"_"$SNa"_out35_TEmatchedCDS.tsv" | sort -t $'\t') <(awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $SNa"_inputGENEs_"$IGEgroup | sort -t $'\t') | awk '!seen[$1]++ {print $5"\t"$6"\t"$7"\t"$1"\t.\t"$8}' > $IGEgroup"_"$SNa"_IGEs.bed"
	# also print a lookup table to link CDS's to the matching TE
	join <(awk '{print $3"\t"$1"\t"$2"\t"$4}' "tmp."$IGEgroup"_"$SNa"_out35_TEmatchedCDS.tsv" | sort -t $'\t') <(awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $SNa"_inputGENEs_"$IGEgroup | sort -t $'\t') | awk '{print $1"\t"$4"\t"$2"\t"$3}' > $IGEgroup"_"$SNa"_IGE_TE_lookup.tsv"
	rm "tmp."$IGEgroup"_"$SNa*
	echo " <-- done finding matching IGE_$IGEgroup at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
}

create_IGE_reference()
{
	echo " --> start creating IGE_$IGEgroup reference files at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
	cd $wd"/IGE_COLLECTION_"$SNa
	mkdir $IGEgroup"_"$SNa"_IGEref_CONTAINER"
	cd $IGEgroup"_"$SNa"_IGEref_CONTAINER"
	# get fasta for CDS's
	bedtools getfasta -fi $REFpath$REFbase".fa" -bed $1 -s -name > $IGEgroup"_"$SNa"_IGEref.TEs.fa"
	# copy original reference genome (this is necessary because RepeatMasker changes current files)
	cp $REFpath$REFbase".fa" "./tmp."$REFbase".fa"
	# run RepeatMasker using newly created CDS fasta file as "TE", then remove unnecessary output
	RepeatMasker -lib $IGEgroup"_"$SNa"_IGEref.TEs.fa" -no_is -nolow -s -pa $nc "./tmp."$REFbase".fa"
	rm "./tmp."$REFbase".fa"
	rm "./tmp."$REFbase".fa.cat.gz"
	rm "./tmp."$REFbase".fa.ori.out"
	rm "./tmp."$REFbase".fa.out"
	rm "./tmp."$REFbase".fa.tbl"
	# output is a "clean" reference genome, now void of IGE sequences
	mv "./tmp."$REFbase".fa.masked" $IGEgroup"_"$SNa"_IGEref.clean.noTEs.fa"
	# create FASTA of IGEs where each IGE has chromosome name ">TEchr_..."
	awk '{if ($1 ~ ">") {gsub(/>/,""); print ">TEchr_"$1"\t"$2"\t"$3} else {print $0}}' $IGEgroup"_"$SNa"_IGEref.TEs.fa" > $IGEgroup"_"$SNa"_IGEref.clean.onlyTEs.fa"
	# combine IGE-cleaned reference genome with IGE-fasta file	
	cat $IGEgroup"_"$SNa"_IGEref.clean.noTEs.fa" $IGEgroup"_"$SNa"_IGEref.clean.onlyTEs.fa" > $IGEgroup"_"$SNa"_IGEref.clean.fa"
	# generate STAR genome index
	mkdir "STAR_"$IGEgroup"_"$SNa"_IGEref"
	STAR --runMode genomeGenerate --genomeFastaFiles $IGEgroup"_"$SNa"_IGEref.clean.fa" --genomeDir "STAR_"$IGEgroup"_"$SNa"_IGEref" --runThreadN $nc
	mv Log.out "STAR_"$IGEgroup"_"$SNa"_IGEref/."
	# generate BLAST databases
	makeblastdb -dbtype nucl -in $IGEgroup"_"$SNa"_IGEref.clean.onlyTEs.fa"
	makeblastdb -dbtype nucl -in $IGEgroup"_"$SNa"_IGEref.clean.noTEs.fa"
	cd $wd
	echo " <-- done creating IGE_$IGEgroup reference files at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
}

align_IGEref_and_filter()
{
	cd $wd"/IGE_COLLECTION_"$SNa
	if [ ! -d $IGEgroup"_"$SNa"_"$SNo"_"$LNo"_IGE" ]; then mkdir $IGEgroup"_"$SNa"_"$SNo"_"$LNo"_IGE"; fi
	cd $IGEgroup"_"$SNa"_"$SNo"_"$LNo"_IGE"
	echo "==================" > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "|| TEchim - IGE || " >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "==================" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " --> start mapping on IGE reference at ... $(date)" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	#### STREAM 1 ####
	# run STAR in chimera-mode
	STAR --runThreadN $nc --genomeDir $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_IGEref_CONTAINER/STAR_"$IGEgroup"_"$SNa"_IGEref" --readFilesIn $1 $2 --chimSegmentMin 20 --chimOutType WithinBAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out4_"
	mkdir $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_STAR"
	mv $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out4"* $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_STAR"/.
	# extract only hits that cross IGE-GENE breakpoints. the awk commands remove IGE-IGE reads
	samtools view $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_STAR"/$IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | grep TEchr_ | awk '($7 != "=") && ($3 !~ "TEchr_" || $7 !~ "TEchr_")' > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5_STREAM1_TExGENES.sam"
	if [ ! -s $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5_STREAM1_TExGENES.sam" ]; then
		echo " #### ERROR: file does not have any TE-GENE brakpoint spanning reads!" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
		echo " --> exited at ... $(date)" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
		exit
	fi
	#### STREAM 2 ####
	# The output of this stream is based on the STAR alignment of in-silico paired-end reads.
	# The hits will only be used for reads where the genome-section is successfully mapped with BLAST (further downstream), but the transposon section is NOT.
	# Using the maximum fragment length, and information about whether the gene- and te- reads are mapped to the (+)ive or (-)ive strand (all contained in SAM-flag) are used for output information.
	# Select pairs where both reads map to positive strand
	samtools view -f 65 -F 48 $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_STAR"/$IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4,$4+mfl,$1"|"$7"|minus|GENE-TE|"s"|"l"|"$8,".","+",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8,$8+mfl,$1"|"$3"|plus|TE-GENE|"s"|"l"|"$4,".","-",$1}}' | sed 's/TEchr_//g' > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5b_STREAM2_FORout10c.bed"
	# Select pairs where read1 maps to positive, and read2 to negative
	samtools view -f 97 -F 16 $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_STAR"/$IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4,$4+mfl,$1"|"$7"|plus|GENE-TE|"s"|"l"|"$8,".","+",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8-mfl,$8,$1"|"$3"|plus|TE-GENE|"s"|"l"|"$4,".","+",$1}}' | sed 's/TEchr_//g' >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5b_STREAM2_FORout10c.bed"
	# Select pairs where read1 maps to negative, and read2 to positive
	samtools view -f 81 -F 32 $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_STAR"/$IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4-mfl,$4,$1"|"$7"|minus|GENE-TE|"s"|"l"|"$8,".","-",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8,$8+mfl,$1"|"$3"|minus|TE-GENE|"s"|"l"|"$4,".","-",$1}}' | sed 's/TEchr_//g' >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5b_STREAM2_FORout10c.bed"
	# Select pairs where both reads map to negative strand
	samtools view -f 113 $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_STAR"/$IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | awk -v s="$SNo" -v l="$LNo"  -v mfl="$MaxFragLength" 'BEGIN {OFS = "\t"} {if ($3 !~ "TEchr_" && $7 ~ "TEchr_") {print $3,$4-mfl,$4,$1"|"$7"|plus|GENE-TE|"s"|"l"|"$8,".","-",$1} else if ($3 ~ "TEchr_" && $7 !~ "TEchr_" && $7 != "=") {print $7,$8-mfl,$8,$1"|"$3"|minus|TE-GENE|"s"|"l"|"$4,".","+",$1}}' | sed 's/TEchr_//g' >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5b_STREAM2_FORout10c.bed"
	bedtools sort -i $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5b_STREAM2_FORout10c.bed" > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out10b_STREAM2_additional_chimera.bed"
	rm $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5b_STREAM2_FORout10c.bed"
	echo " ------ sample contains $(wc -l $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam" | awk '{print $1}') reads that span gene-TE breakpoint." >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " <-- done with mapping on IGE reference at ... $(date)" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "--------------------------------" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	cd $wd
}

blast_on_longreads ()
{
	cd $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_IGE"
	echo " --> start BLAST at ... $(date)" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# extract readnames and remove duplicates
	awk '{print $1}' $1 | sort | uniq > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt"
	# combine readnames with long sequences stored in the lookup file
	# LANG=en_EN is a bug-fix to make sort compatible with join
	LANG=en_EN join -1 1 -2 1 <(sort $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt") <(zcat < $2 ) > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv" && rm $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt"
	# the following while loop will add data to file using ">>". just as safety
	# precaution, this line makes sure no data with such a name exists
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
	# loop through long sequences. two blast searches are performed, one on a
	# reference genome without TE sequences and another one on all TE sequences.
	# only the top hit is reported.
	while IFS=$' ' read -r readname sequence
	do	
		echo "${readname}" > var1
		echo "${sequence}" > var2
		echo "${sequence}" | blastn -db $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_IGEref_CONTAINER/"$IGEgroup"_"$SNa"_IGEref.clean.noTEs.fa" -outfmt "6 qstart qend sseqid sstart send sstrand qlen" -num_alignments 1 -num_threads $nc | head -n 1 > var3		
		if ! [ -s var3 ]; then echo "noGENEfound" > var3; fi
		echo "${sequence}" | blastn -db $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_IGEref_CONTAINER/"$IGEgroup"_"$SNa"_IGEref.clean.onlyTEs.fa" -outfmt "6 qstart qend sseqid sstart send slen sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var4
		if ! [ -s var4 ]; then echo "noTEfound" > var4; fi
		paste var1 var2 var3 var4 >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
		rm var*
	done < $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv"
	# remove reads that did not give BLAST result for genomic location
	grep -v noGENEfound $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv" > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv"
	# extract reads for which TE has been identified
	grep -v noTEfound $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv" | awk 'BEGIN {OFS = "\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16}' > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out10a_STREAM1_TExGENES_blastedreads.tsv"
	# extract reads where NO TE was found - these are filtered (at least $fastalength/2 should NOT map to genome) and STREAM2 information is added
	# remove reads where less than half the $fastalength is left unmapped when BLASTed to no_TE_genome	
	grep noTEfound $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv" | awk -v flength="$fastalength" '{if ($9-$4 > flength/2 || $3 > flength/2) {print $0} }' > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9b_fromSTREAM1forSTREAM2_TExGENES_findTE.tsv"
	join -1 1 -2 7 <(sort $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9b_fromSTREAM1forSTREAM2_TExGENES_findTE.tsv") <(sort -k 7 $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out10b_STREAM2_additional_chimera.bed") | awk 'BEGIN {OFS = "\t"} {if ($8=="plus") {if ($14 ~ "[\|]TE-GENE[\|]") {print $5,$6-1,$6,$14,".","+",$3,$4,$9} else if ($14 ~ "[\|]GENE-TE[\|]") {print $5,$7-1,$7,$14,".","+",$3,$4,$9}} else if ($8=="minus") {if ($14 ~ "[\|]TE-GENE[\|]") {print $5,$6-1,$6,$14,".","-",$3,$4,$9} else if ($14 ~ "[\|]GENE-TE[\|]") {print $5,$7-1,$7,$14,".","-",$3,$4,$9}}}' | bedtools sort -i - > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11b_STREAM2_additional_chimera.filtered.bed"
	tr "|" "\t" < $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11b_STREAM2_additional_chimera.filtered.bed" | awk -v flength="$fastalength" 'BEGIN {OFS="\t"} {if ($6 == "plus") {if ($7 == "TE-GENE") {testart=$10+flength; teend=$10+$13} else if ($7 == "GENE-TE") {testart=$10+flength+$14-$15; teend=$10}} else if ($6 == "minus") {if ($7 == "TE-GENE") {testart=$10+flength-$13; teend=$10} else if ($7 == "GENE-TE") {testart=$10+flength; teend=$10+$15-$14}} {print $1,$2,$3,$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"testart"-"teend"|0",$11,$12}}' > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed"
	echo " ------ BLAST results:" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " ------ $(wc -l $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv" | awk '{print $1}') input reads" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " ------ $(wc -l $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv" | awk '{print $1}') reads with a genome location" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " ------ |> of the reads with genome location: $(wc -l $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out10a_STREAM1_TExGENES_blastedreads.tsv" | awk '{print $1}') reads with a genome- and transposon location" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " ------ |> of the reads with genome location: $(wc -l $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed" | awk '{print $1}') reads from STREAM2 will be added to STREAM1" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"	
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv"
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads_plusTEnohits.tsv"
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out9b_fromSTREAM1forSTREAM2_TExGENES_findTE.tsv"
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out10b_STREAM2_additional_chimera.bed"
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11b_STREAM2_additional_chimera.filtered.bed"
	echo " <-- done with BLAST at ... $(date)" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "--------------------------------" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	cd $wd
}

create_summary_table ()
{
	cd $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_IGE"
	echo " --> start creating summary table at ... $(date)" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# determine whether the section that maps to the genome is the 5' or the 3' end
	# of the mRNA section:
	# 5'-######|GENE|TE|######-3' => PART5
	# 5'-######|TE|GENE|######-3' => PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; if (a < b) {print "PART5"} else {print "PART3"}}' < $1 > tmpfile.genepart
	# determine the precise breakpoint on the chromosome. this depends on  whether
	# the chromosomal part is PART5 or PART3
	awk 'BEGIN {OFS = "\t"} {a = $5; {print a}}' < $1 > tmpfile.chr
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d-1} else {print c-1}}' < $1 > tmpfile.breakpoint.chr.start
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d} else {print c}}' < $1 > tmpfile.breakpoint.chr.end
	# determine the precise breakpoint on the TE. this depends on  whether the TE
	# part is PART5 or PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $13 ; d = $12 ; if (a < b) {print d} else {print c}}' < $1 > tmpfile.breakpoint.TE
	# determine the overlap between the two mapped sections of the long read
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $4 ; d = $10 ; if (a < b) {print b-c-1} else {print a-d-1}}' < $1 > tmpfile.uncertainty
	awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
		a = $11
		gsub(/TEchr_/,"",a)
		b = $15
		c = $3
		d = $9
		e = $1
		if (c < d) {print e"|"a"|"b"|GENE-TE|"s"|"l} else {print e"|"a"|"b"|TE-GENE|"s"|"l}
		}' < $1 > tmpfile.readname
	awk 'BEGIN {OFS = "\t"} {a = $8 ; if (a == "plus") {print "+"} else {print "-"}}' < $1 > tmpfile.breakpoint.chr.strand
	awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "."}}' < $1 > tmpfile.score
	paste -d'|' tmpfile.readname tmpfile.breakpoint.TE tmpfile.uncertainty > tmpfile.readname.extended
	paste -d'\t' tmpfile.chr tmpfile.breakpoint.chr.start tmpfile.breakpoint.chr.end tmpfile.readname.extended tmpfile.score tmpfile.breakpoint.chr.strand > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11a_STREAM1_breakpoints.bed"
	# add results from STREAM 2
	cat $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11a_STREAM1_breakpoints.bed" $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed" | bedtools sort -i - > $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out13_breakpoints.bed"
	echo " ------ In total, $(wc -l $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out13_breakpoints.bed" | awk '{print $1}') chimeric reads were found." >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm tmpfile.*
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11a_STREAM1_breakpoints.bed"
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out11b1_STREAM2_additional_chimera.filtered.bed"
	rm -f $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out12_STREAM1and2_breakpoints.bed"
	echo " <-- all done at ... $(date)" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "================================" >> $IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm $1
	cd $wd
}

process_P1out_IGE()
{
	cd $wd"/IGE_COLLECTION_"$SNa
	echo " --> start processing PART1 output for IGE_$IGEgroup at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"	
	cat $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_"*"_IGE/"$IGEgroup"_"$SNa"_IGEref_"*"_out13_breakpoints.bed" | bedtools sort -i - > $IGEgroup"_"$SNa"_IGEref_in10_combined.sorted.bed"
	if [ $stranded = "0" ]; then
		bedtools intersect -wa -a $IGEgroup"_"$SNa"_IGEref_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj > $IGEgroup"_"$SNa"_IGEref_out01_genetagged.tsv"
	else
		bedtools intersect -wa -a $IGEgroup"_"$SNa"_IGEref_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj -s > $IGEgroup"_"$SNa"_IGEref_out01_genetagged.tsv"
	fi
	# filter reads inside same gene
	awk '{if ($4 !~ $10) print $0}' $IGEgroup"_"$SNa"_IGEref_out01_genetagged.tsv" > $IGEgroup"_"$SNa"_IGEref_out01a_filtered_genetagged.tsv"
	# Also extract those reads that were not inside annotated gene
	#grep -Fvf <(tr "|" "\t" < $IGEgroup"_"$SNa"_IGEref_in10_combined.sorted.bed" | awk '{print $4}') <(cat $IGEgroup"_"$SNa"_IGEref_out01_genetagged.tsv") > $IGEgroup"_"$SNa"_IGEref_out01b_ReadsOutsideGenes.tsv"
	# separate | delimited field
	tr '|' '\t' < $IGEgroup"_"$SNa"_IGEref_out01a_filtered_genetagged.tsv" > $IGEgroup"_"$SNa"_IGEref_out02_sepparated.tsv" && rm $IGEgroup"_"$SNa"_IGEref_out01a_filtered_genetagged.tsv" && rm $IGEgroup"_"$SNa"_IGEref_out01_genetagged.tsv"
	# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
	# and create column with BASE readname (no :A or :B)
	awk 'BEGIN {OFS = "\t"} {
		a = $5
		gsub(/_LTR/,"",a)
		c = substr($4, 1, length($4)-2)
		print $0"\t"a"\t"c
		}' < $IGEgroup"_"$SNa"_IGEref_out02_sepparated.tsv" > $IGEgroup"_"$SNa"_IGEref_out03pre_TEbase.tsv" && rm $IGEgroup"_"$SNa"_IGEref_out02_sepparated.tsv"
	wc_predup=$(wc -l $IGEgroup"_"$SNa"_IGEref_out03pre_TEbase.tsv" | awk '{print $1}')
	# remove duplicate reads (where both :A and :B version of the same reads were picked up
	awk '!seen[$21]++' $IGEgroup"_"$SNa"_IGEref_out03pre_TEbase.tsv" > $IGEgroup"_"$SNa"_IGEref_out03_TEbase.tsv" && rm $IGEgroup"_"$SNa"_IGEref_out03pre_TEbase.tsv"
	echo " ---- $(wc -l $IGEgroup"_"$SNa"_IGEref_out03_TEbase.tsv" | awk '{print $1}') reads ($wc_predup before removing duplicates)" >> $wd"/"$SNa"_IGE_"$logname".log"
	echo " <-- done processing PART1 output for IGE_$IGEgroup at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
}

combine_hits_of_each_IGE()
{
	cd $wd"/IGE_COLLECTION_"$SNa
	echo " --> start combining IGE_$IGEgroup reads at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
	# create list of all TEs in dataset
	cut -f20 $1 | sort | uniq > $IGEgroup"_"$SNa"_IGEref_out03a_uniqueTEs.tsv"
	# to create collection file, first make sure this file does not yet exist
	rm -f $IGEgroup"_"$SNa"_IGEref_out12_OUTPUT.tsv"
	# loop through all TEs in the dataset
	while read TE
	do
		# grep TEs from combined data set
		awk -v t="$TE" '{if ($20 == t) print $0;}' $1 > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out04.tsv"
		if [ -s "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out04.tsv" ]; then	
			# SPLIT into TE-GENE and GENE-TE
			awk '{if ($7 == "TE-GENE") print $0;}' "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out04.tsv" > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv"
			# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
			# header: Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllBreakpoints(TE)
			if [ -s "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv" ]; then
				bedtools merge -i "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,6,7,17,8,10,5,21 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse,collapse,distinct,collapse -d 20 > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out06_TE-GENE.tsv"
			fi
			awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out04.tsv" > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out05_GENE-TE.tsv"
			if [ -s "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out05_GENE-TE.tsv" ]; then
				bedtools merge -i "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,6,7,17,8,10,5,21 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse,collapse,distinct,collapse -d 20 > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out06_GENE-TE.tsv"
			fi
			# COMBINE TE-GENE and GENE-TE
			cat "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out06_"*".tsv" > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out08.tsv"
			bedtools sort -i "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out08.tsv" > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out09.tsv"
			# create list of all genes
			# first, replace "," with tab, then sort, then pick unique strings, then remove ""
			cut -f12 "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out09.tsv" | tr ',' '\n' | sort | uniq | sed '/\./d' > "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out11_genes.tsv"
			# this list of genes is used to add non-anchor hits to ANCHOR hits, if they overlap with the same gene and TE
			while read GENE
			do
				# create output file that contains for every gene all the insertions of that TE
				if [ $stranded = "0" ]; then
					awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$7"\t"$8"\tunstranded\t"$13"\t"$14"\t"$15"\t"$16;}' "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out09.tsv" >> $IGEgroup"_"$SNa"_IGEref_out12_OUTPUT.tsv"
				else
					awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$7"\t"$8"\tmRNA\t"$13"\t"$14"\t"$15"\t"$16;}' "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out09.tsv" >> $IGEgroup"_"$SNa"_IGEref_out12_OUTPUT.tsv"
				fi			
			done < "tmp."$IGEgroup"_"$SNa"_"$TE"_IGEref_out11_genes.tsv"
			rm -f "tmp."$IGEgroup"_"$SNa"_"$TE*
		fi
	done < $IGEgroup"_"$SNa"_IGEref_out03a_uniqueTEs.tsv"
	# remove empty lines
	grep . $IGEgroup"_"$SNa"_IGEref_out12_OUTPUT.tsv" > $IGEgroup"_"$SNa"_IGEref_out12a_OUTPUT.tsv"
	rm $IGEgroup"_"$SNa"_IGEref_out03a_uniqueTEs.tsv"
	rm $IGEgroup"_"$SNa"_IGEref_out03_TEbase.tsv"
	rm $IGEgroup"_"$SNa"_IGEref_out12_OUTPUT.tsv"
	echo " <-- done combining IGE_$IGEgroup reads at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
}

check_for_IGE_splicesites()
{
	cd $wd"/IGE_COLLECTION_"$SNa
	echo " --> start checking for IGE_$IGEgroup breaks at splice sites at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
	# now, every line is one incident of a gene-TE chimera
	while read line
	do
		# convert line to bed format
		# header: Chr(genome)\tStart(genome)[this is the breakpoint minus 1]\tEnd(genome)[this is the breakpoint]\tTE|Gene|TE-GENEorGENE-TE\t"."\tStrand(genome)
		echo $line | awk '{ print $2"\t"$4-1"\t"$4"\t"$5"|"$1"|"$6"\t.\t"$3}' > "tmp."$IGEgroup"_"$SNa"_IGEref_out13.bed"
		# $a is TE-GENE or GENE-TE	
		a=$(echo $line | awk '{print $6}')
		# $b is the gene
		b=$(echo $line | awk '{print $1}')
		# intersect with FEATURES file to get all ovelapping features (output is ;-separated)
		if [ $stranded = "0" ]; then
			d=$(bedtools intersect -a "tmp."$IGEgroup"_"$SNa"_IGEref_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa | awk '{print $10}' | paste -sd ";" -)
		else
			d=$(bedtools intersect -a "tmp."$IGEgroup"_"$SNa"_IGEref_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa | awk '{print $10}' | paste -sd ";" -)
		fi
		# check whether breakpoint on the genome overlaps with splice donor site (in the case of a GENE-TE fragment)
		# or with splice acceptor site (in the case of a TE-GENE fragment). the REF files _SPLICE_DONORS.bed and 
		# _SPLICE_ACCEPTORS.bed have a +/- 10nt, but the precise breakpoint is extracted here (which is the sequence
		# number of the last nucleotide that is still part of the exon for SPLICE_DONORS and the first exonic nt 
		# for SPLICE_ACCEPTORS.
		# for GENE-TE fragments, the exon is the splice donor
		if [ $a == "GENE-TE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_DONORS.bed" > "tmp."$IGEgroup"_"$SNa"_IGEref_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$IGEgroup"_"$SNa"_IGEref_out13.bed" -b "tmp."$IGEgroup"_"$SNa"_IGEref_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$IGEgroup"_"$SNa"_IGEref_out14_ref_for_overlap.bed"
		# for TE-GENE fragments, the exon is the splice acceptor
		elif [ $a == "TE-GENE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_ACCEPTORS.bed" > "tmp."$IGEgroup"_"$SNa"_IGEref_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$IGEgroup"_"$SNa"_IGEref_out13.bed" -b "tmp."$IGEgroup"_"$SNa"_IGEref_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$IGEgroup"_"$SNa"_IGEref_out14_ref_for_overlap.bed"
		else
			c="N/A"
		fi
		if [[ $c = "" ]]; then c="."; fi
		# the loop aboce created the unix variable $c. here, this variable is appended to the main file.
		echo $line | awk -v c="$c" -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"c"\t"d}'
		rm "tmp."$IGEgroup"_"$SNa"_IGEref_out13.bed"
	done < $1 > $IGEgroup"_"$SNa"_IGEref_out15.tsv"
	rm $1
	echo " <-- done checking for IGE_$IGEgroup breaks at splice sites at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
}

split_IGE_breakpoints()
{
	cd $wd"/IGE_COLLECTION_"$SNa
	echo " --> start tyding up IGE_$IGEgroup breakpoints at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
	# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)
	cat $1 | while read line; do echo $line | awk '{print $12}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| awk '{print $0"\n"}' | awk 'NF'; done > $IGEgroup"_"$SNa"_IGE_newcolb"
	cat $1 | while read line; do echo $line | awk '{print $11}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| awk '{print $0"\n"}' | awk 'NF'; done > $IGEgroup"_"$SNa"_IGE_newcolc"
	# for final output, get rid of original TE-breakpoint field
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$13,$14,$15,$16}' $1 > $IGEgroup"_"$SNa"_IGE_newcola"
	# append pooled TE breakpoints to final output
	paste $IGEgroup"_"$SNa"_IGE_newcola" $IGEgroup"_"$SNa"_IGE_newcolb" $IGEgroup"_"$SNa"_IGE_newcolc" > $IGEgroup"_"$SNa"_IGEref_chimericreads_final.tsv"
	awk '{all=$0; gsub(/@/,"\t"); if($5 !~ $16) print all}' $IGEgroup"_"$SNa"_IGEref_chimericreads_final.tsv" > $IGEgroup"_"$SNa"_IGEref_chimericreads_final.FILTERED.tsv"
	rm $IGEgroup"_"$SNa"_IGE_newcola" $IGEgroup"_"$SNa"_IGE_newcolb"
	rm $1
	echo " <-- done tyding up IGE_$IGEgroup breakpoints at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
}


################################################################################
################################################################################
# execute:

# change to wd
cd $wd
																																		     
get_variables
write_logfile
split_CDS
calculate_TE_coverage
input_letters=$(printf "a\tb\tc\td\te\tf\tg\th\ti\tj")	
for IGEgroup in $input_letters
do
	find_matching_IGEs
	create_IGE_reference $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_IGEs.bed"
	list_of_SNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}' | rev)
	for SNo in $list_of_SNo
	do
		list_of_LNo=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | awk '{gsub(/_/,"\t"); print $2"\t"$1}' | rev | grep $SNo | awk '{print $1}')
		for LNo in $list_of_LNo
		do
			(if [[ ! -f $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_LOOKUP.fa" ]]; then zcat < $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_LOOKUP.sorted.tsv.gz" | awk '{print ">"$1"\n"$2}' > $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_LOOKUP.fa"; fi
			align_IGEref_and_filter $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_out3_1.trimmed.fq" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_out3_2.trimmed.fq" 
			blast_on_longreads $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_IGE/"$IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out5_STREAM1_TExGENES.sam" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_LOOKUP.sorted.tsv.gz"
			create_summary_table $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_"$SNo"_"$LNo"_IGE/"$IGEgroup"_"$SNa"_IGEref_"$SNo"_"$LNo$"_out10a_STREAM1_TExGENES_blastedreads.tsv") &
		done
	done
	wait
	process_P1out_IGE		
	combine_hits_of_each_IGE $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_IGEref_out03_TEbase.tsv"
	check_for_IGE_splicesites $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_IGEref_out12a_OUTPUT.tsv"
	split_IGE_breakpoints $wd"/IGE_COLLECTION_"$SNa"/"$IGEgroup"_"$SNa"_IGEref_out15.tsv"
done
echo " <-- all done at ... $(date)" >> $wd"/"$SNa"_IGE_"$logname".log"
