#!/bin/bash

################################################################################
# TITLE: TEchim - PART 3
# VERSION: 0.2.0 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 04/06/2019 (dd/mm/yyyy)
# DESCRIPTION: This script takes the output from PART1 to 1. calculate the
# average TE expression levels 2. to find matching immobile CDS's with similar
# expression levels and 3. to use these CDS's as IGEs to generate a sample-
# specific IGE reference genome. 
################################################################################

################################################################################
# REQUIREMENTS:
# - STAR (https://github.com/alexdobin/STAR)
# - samtools
# - bedtools
################################################################################

################################################################################
################################################################################
# set parameters
wd=~							# working directory
path_to_PART1_output=/PATH/TO/OUTPUT/
SNa=NAME_OF_EXP					# sample name	
nc=1							# number of cores
REFpath=/PATH/TO/REF/			# same path as in _buildREF
REFbase=dmel625					# same name as in _buildREF
################################################################################
################################################################################
# functions:

find_matching_IGEs()
{
	echo " --> start finding matching IGEs at ... $(date)" >> $SNa"_PART3_"$logname".log"
	# randomly subsample filtered CDS (here: take 10%)
	tenpc_CdS=$(wc -l $REFpath$REFbase".CDS_for_IGE.bed" | awk '{OFMT = "%.0f"; print $1/10}')
	sort -R $REFpath$REFbase".CDS_for_IGE.bed" | head -n $tenpc_CdS > $SNa"_inputGENEs.bed"
	# assess how many Samples and Lanes exist - run bedtools mutlicov on each sample/lane
	list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}')
	for snum in $list_of_snum
	do
		list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2"\t"$3}' | grep $snum | awk '{print $2}')
		for lnum in $list_of_lanes
		do
			# check if .bam file is indexed, if not then index it
			if [[ ! -f $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam.bai" ]]; then samtools index $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam"; fi
			# quantify coverage of TEs in each sample/lane - this line will print the number
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $REFpath"TEs.fa.bed" -s | awk '{print $7}' > "tmp."$SNa"_"$snum"_"$lnum"_out31_TEcoverage_perSample.tsv"
			# this line will print the TE name - only has to be done once (hence the if statement)
			if [[ ! -f "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv" ]]; then bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $REFpath"TEs.fa.bed" -s | awk '{print $4}' > "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv"; fi
			# quantify coverage of CDS (subsample) in each sample/lane - this line will print the number
			bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $SNa"_inputGENEs.bed" -s | awk '{print $7}' > "tmp."$SNa"_"$snum"_"$lnum"_out31_GENEcoverage_perSample.tsv"
			# this line will print the GENE names - only has to be done once (hence the if statement)
			if [[ ! -f "tmp."$SNa"_out31a_GENEcoverage_GENEnames.tsv" ]]; then bedtools multicov -bams $path_to_PART1_output$SNa"_"$snum"_"$lnum"/"$SNa"_"$snum"_"$lnum"_STAR/"$SNa"_"$snum"_"$lnum"_out4_Aligned.sortedByCoord.out.bam" -bed $SNa"_inputGENEs.bed" -s | awk '{print $4}' > "tmp."$SNa"_out31a_GENEcoverage_GENEnames.tsv"; fi	
		done
	done
	# combine all TE coverage numbers					  | calculate average, round to full integer
	paste "tmp."$SNa"_"*"_out31_TEcoverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$SNa"_out31b_TEcoverage_averages.tsv"
	# combine averages with TE names
	paste "tmp."$SNa"_out31a_TEcoverage_TEnames.tsv" "tmp."$SNa"_out31b_TEcoverage_averages.tsv" | sort -k2n > "tmp."$SNa"_out32_TEcoverage_averages.tsv"
	# combine all GENE coverage numbers						| calculate average, round to full integer
	paste "tmp."$SNa"_"*"_out31_GENEcoverage_perSample.tsv" | awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; printf "%.0f\n",sum}' > "tmp."$SNa"_out31b_GENEcoverage_averages.tsv"
	# combine averages with GENE names
	paste "tmp."$SNa"_out31a_GENEcoverage_GENEnames.tsv" "tmp."$SNa"_out31b_GENEcoverage_averages.tsv" | sort -k2n > "tmp."$SNa"_out32_GENEcoverage_averages.tsv"
	# loop through all TEs
	while read TEline
	do
		coverage=$(echo $TEline | awk '{print $2}')
		# find the gene that matches the expression level of each TE									  | print matching TE-GENE pairs
		awk -v v1="$coverage" '$2>v1 {print $1"\t"$2; exit}' "tmp."$SNa"_out32_GENEcoverage_averages.tsv" | awk -v te="$TEline" '{print te"\t"$0}'
	done < "tmp."$SNa"_out32_TEcoverage_averages.tsv" > "tmp."$SNa"_out35_TEmatchedCDS.tsv"
	# left-join TE-GENE pairs with original GENE bed (this is to get the strand of the according CDS's)																				 | print a bed file with each UNIQUE CDS (strand-specific!)
	join <(awk '{print $3"\t"$1"\t"$2"\t"$4}' "tmp."$SNa"_out35_TEmatchedCDS.tsv" | sort -t $'\t') <(awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $SNa"_inputGENEs.bed" | sort -t $'\t') | awk '!seen[$1]++ {print $5"\t"$6"\t"$7"\t"$1"\t.\t"$8}' > $SNa"_IGEs.bed"
	# also print a lookup table to link CDS's to the matching TE
	join <(awk '{print $3"\t"$1"\t"$2"\t"$4}' "tmp."$SNa"_out35_TEmatchedCDS.tsv" | sort -t $'\t') <(awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' $SNa"_inputGENEs.bed" | sort -t $'\t') | awk '{print $1"\t"$4"\t"$2"\t"$3}' > $SNa"_IGE_TE_lookup.tsv"
	rm "tmp."$SNa*
	echo " <-- done finding matching IGEs at ... $(date)" >> $SNa"_PART3_"$logname".log"
}

create_IGE_reference()
{
	echo " --> start creating IGE reference files at ... $(date)" >> $SNa"_PART3_"$logname".log"
	# get fasta for CDS's
	bedtools getfasta -fi $REFpath$REFbase".fa" -bed $1 -s -name > $SNa"_IGEref.TEs.fa"
	# copy original reference genome (this is necessary because RepeatMasker changes current files)
	cp $REFpath$REFbase".fa" "./tmp."$REFbase".fa"
	# run RepeatMasker using newly created CDS fasta file as "TE", then remove unnecessary output
	RepeatMasker -lib $SNa"_IGEref.TEs.fa" -no_is -nolow -s -pa $nc "./tmp."$REFbase".fa"
	rm "./tmp."$REFbase".fa"
	rm $REFbase".fa.cat.gz"
	rm $REFbase".fa.ori.out"
	rm $REFbase".fa.out"
	rm $REFbase".fa.tbl"
	# output is a "clean" reference genome, now void of IGE sequences
	mv $REFbase".fa.masked" $SNa"_IGEref.clean.noTEs.fa"
	# create FASTA of IGEs where each IGE has chromosome name ">TEchr_..."
	awk '{if ($1 ~ ">") {gsub(/>/,""); print ">TEchr_"$1"\t"$2"\t"$3} else {print $0}}' $SNa"_IGEref.TEs.fa" > $SNa"_IGEref.clean.onlyTEs.fa"
	# combine IGE-cleaned reference genome with IGE-fasta file	
	cat $SNa"_IGEref.clean.noTEs.fa" $SNa"_IGEref.clean.onlyTEs.fa" > $SNa"_IGEref.clean.fa"
	# generate STAR genome index
	mkdir "STAR_"$SNa"_IGEref"
	STAR --runMode genomeGenerate --genomeFastaFiles $SNa"_IGEref.clean.fa" --genomeDir "STAR_"$SNa"_IGEref" --runThreadN $nc
	mv Log.out "STAR_"$SNa"_IGEref/."
	# generate BLAST databases
	makeblastdb -dbtype nucl -in $SNa"_IGEref.clean.onlyTEs.fa"
	makeblastdb -dbtype nucl -in $SNa"_IGEref.clean.noTEs.fa"
	echo " <-- done creating IGE reference files at ... $(date)" >> $SNa"_PART3_"$logname".log"
}

################################################################################
################################################################################
# execute:

# change to wd
cd $wd
mkdir $SNa"_IGEref_CONTAINER"
cd $SNa"_IGEref_CONTAINER"

# create .log file
logname=$(date | awk '{gsub(/\:/,"-",$5); print $4$3$2"_"$5}')
echo "====================" > $SNa"_PART3_"$logname".log"
echo "|| TEchim - PART3 || " >> $SNa"_PART3_"$logname".log"
echo "====================" >> $SNa"_PART3_"$logname".log"
echo "Parameters:" >> $SNa"_PART3_"$logname".log"
echo "Working directory:" "$wd/$SNa"_IGEref_CONTAINER"" >> $SNa"_"$logname".log"
echo "Location of PART1 output:" "$path_to_PART1_output" >> $SNa"_"$logname".log"
echo "Sample name:" "$SNa" >> $SNa"_PART3_"$logname".log"
echo "Reference files:" "$REFpath$REFbase" >> $SNa"_PART3_"$logname".log"
echo "--------------------------------" >> $SNa"_PART3_"$logname".log"
echo " --> starting at ... $(date)" >> $SNa"_PART3_"$logname".log"
echo "--------------------------------" >> $SNa"_PART3_"$logname".log"

find_matching_IGEs

create_IGE_reference $SNa"_IGEs.bed"
