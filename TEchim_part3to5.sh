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
wd=~							# working directory (no trailing "/")
path_to_PART1_output=/PATH/TO/OUTPUT/	# with trailing "/"
SNa=NAME_OF_EXP					# sample name	
nc=1							# number of cores
REFpath=/PATH/TO/REF/			# same path as in _buildREF (with trailing "/")
REFbase=dmel625					# same name as in _buildREF
################################################################################
################################################################################
# functions:

find_matching_IGEs()
{
	echo " --> start finding matching IGEs at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
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
	echo " <-- done finding matching IGEs at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
}

create_IGE_reference()
{
	echo " --> start creating IGE reference files at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
	mkdir $wd"/"$SNa"_IGEref_CONTAINER"
	cd $wd"/"$SNa"_IGEref_CONTAINER"
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
	cd $wd
	echo " <-- done creating IGE reference files at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
}

align_and_filter()
{
	if [ ! -d $wd"/IGE_"$SNa"_"$SNo"_"$LNo ]; then mkdir $wd"/IGE_"$SNa"_"$SNo"_"$LNo; fi
	cd $wd"/IGE_"$SNa"_"$SNo"_"$LNo
	echo " --> start mapping at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# run STAR in chimera-mode
	STAR --runThreadN $nc \
		--genomeDir $wd"/"$SNa"_IGEref_CONTAINER/STAR_"$SNa"_IGEref" \
		--readFilesIn $1 $2 \
		--chimSegmentMin 20 \
		--chimOutType WithinBAM \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix $SNa"_IGEref_"$SNo"_"$LNo$"_out4_"	
	# extract only hits that cross IGE-GENE breakpoints. the awk commands remove IGE-IGE reads
	samtools view $SNa"_IGEref_"$SNo"_"$LNo$"_out4_"Aligned.sortedByCoord.out.bam | grep TEchr_ | awk '($7 != "=")' > $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam"
	mkdir $SNa"_IGEref_"$SNo"_"$LNo"_STAR"
	mv $SNa"_IGEref_"$SNo"_"$LNo$"_out4"* $SNa"_IGEref_"$SNo"_"$LNo"_STAR"/.
	if [ ! -s $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam" ]; then
		echo " #### ERROR: file does not have any TE-GENE brakpoint spanning reads!" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
		echo " --> exited at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
		exit
	fi
	echo " ------ sample contains $(wc -l $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam" | awk '{print $1}') reads that span gene-TE breakpoint." >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo " <-- done with mapping at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "--------------------------------" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	cd $wd
}

blast_on_longreads ()
{
	cd $wd"/IGE_"$SNa"_"$SNo"_"$LNo
	echo " --> start BLAST at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# extract readnames and remove duplicates
	awk '{print $1}' $1 | sort | uniq > $SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt"
	# combine readnames with long sequences stored in the lookup file
	# LANG=en_EN is a bug-fix to make sort compatible with join
	LANG=en_EN join -1 1 -2 1 <(sort $SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt") <(zcat < $2 ) > $SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv" && rm $SNa"_IGEref_"$SNo"_"$LNo$"_out6_TExGENES_readnames.txt"
	# the following while loop will add data to file using ">>". just as safety
	# precaution, this line makes sure no data with such a name exists
	rm -f $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
	# loop through long sequences. two blast searches are performed, one on a
	# reference genome without TE sequences and another one on all TE sequences.
	# only the top hit is reported.
	while IFS=$' ' read -r readname sequence
	do	
		echo "${readname}" > var1
		echo "${sequence}" > var2
		echo "${sequence}" | blastn -db $wd"/"$SNa"_IGEref_CONTAINER/"$SNa"_IGEref.clean.noTEs.fa" -outfmt "6 qstart qend sseqid sstart send sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var3		
		if ! [ -s var3 ]; then echo "no-hit" > var3; fi
		echo "${sequence}" | blastn -db $wd"/"$SNa"_IGEref_CONTAINER/"$SNa"_IGEref.clean.onlyTEs.fa" -outfmt "6 qstart qend sseqid sstart send slen sstrand" -num_alignments 1 -num_threads $nc | head -n 1 > var4
		if ! [ -s var4 ]; then echo "no-hit" > var4; fi
		
		##### HERE: if readname sequence == blastn sequence, then discard
		
		paste var1 var2 var3 var4 >> $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
		rm var*
	done < $SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv"
	# remove reads that did not give BLAST result
	grep -v no-hit $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv" > $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv"
	echo " ------ BLAST results: $(wc -l $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" | awk '{print $1}') hits ($(grep no-hit $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv" | wc -l | awk '{print $1}') no hit)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm -f $SNa"_IGEref_"$SNo"_"$LNo$"_out7_TExGENES_longreads.tsv"
	rm -f $SNa"_IGEref_"$SNo"_"$LNo$"_out8_TExGENES_blastedreads_plusnohit.tsv"
	echo " <-- done with BLAST at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "--------------------------------" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam"
	cd $wd
}

create_summary_table ()
{
	cd $wd"/IGE_"$SNa"_"$SNo"_"$LNo
	echo " --> start creating summary table at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	# determine whether the section that maps to the genome is the 5' or the 3' end
	# of the mRNA section:
	# 5'-######|GENE|TE|######-3' => PART5
	# 5'-######|TE|GENE|######-3' => PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; if (a < b) {print "PART5"} else {print "PART3"}}' < $1 > tmpfile.genepart
	# determine the precise breakpoint on the chromosome. this depends on  whether
	# the chromosomal part is PART5 or PART3
	awk 'BEGIN {OFS = "\t"} {a = $5; {print a}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.chr
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d-1} else {print c-1}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.start
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $6 ; d = $7 ; if (a < b) {print d} else {print c}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.end
	# determine the precise breakpoint on the TE. this depends on  whether the TE
	# part is PART5 or PART3
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $13 ; d = $12 ; if (a < b) {print d} else {print c}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.TE
	# determine the overlap between the two mapped sections of the long read
	awk 'BEGIN {OFS = "\t"} {a = $3 ; b = $9 ; c = $4 ; d = $10 ; if (a < b) {print b-c-1} else {print a-d-1}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.uncertainty
	if [[ $stranded = "0" ]]; then
		awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
			a = $11
			gsub(/TEchr_/,"",a)
			if ($8 == "plus") {if ($15=="plus") {b = "forward"} else {b = "reverse"}} else {if ($15 == "plus") { b = "reverse" } else { b = "forward" }}
			if ($3 < $9) {print $1"|"a"|"b"|GENE-TE|S"s"|L"l} else {print $1"|"a"|"b"|TE-GENE|S"s"|L"l}
			}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.readname
			awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "+"}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.strand
	else
	awk -v s="$SNo" -v l="$LNo" 'BEGIN {OFS = "\t"} {
		a = $11
		gsub(/TEchr_/,"",a)
		b = $15
		c = $3
		d = $9
		e = $1
		if (c < d) {print e"|"a"|"b"|GENE-TE|S"s"|L"l} else {print e"|"a"|"b"|TE-GENE|S"s"|L"l}
		}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.readname
		awk 'BEGIN {OFS = "\t"} {a = $8 ; if (a == "plus") {print "+"} else {print "-"}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.breakpoint.chr.strand
	fi
	awk 'BEGIN {OFS = "\t"} {a = $1 ; {print "."}}' < $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" > tmpfile.score
	paste -d'|' tmpfile.readname tmpfile.breakpoint.TE tmpfile.uncertainty > tmpfile.readname.extended
	paste -d'\t' tmpfile.chr tmpfile.breakpoint.chr.start tmpfile.breakpoint.chr.end tmpfile.readname.extended tmpfile.score tmpfile.breakpoint.chr.strand > $SNa"_IGEref_"$SNo"_"$LNo$"_out10_breakpoints.bed"
	#paste -d'\t' $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" tmpfile.breakpoint.chr.end tmpfile.breakpoint.TE tmpfile.uncertainty > $SNa"_IGEref_"$SNo"_"$LNo$"_out11_combined_results.tsv" &&
	rm tmpfile.*
	echo " <-- all done at ... $(date)" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	echo "================================" >> $SNa"_IGEref_"$SNo"_"$LNo"_"$logname".log"
	rm $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv"
	cd $wd
}

process_PART1_output()
{
	cd $wd
	echo " --> start processing PART1 output at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
	# filter reads inside same gene
	awk '{if ($4 !~ $10) print $0}' $1 > $SNa"_IGEref_out01a_filtered_genetagged.tsv"
	# separate | delimited field
	sed -e $'s/|/\t/g' $SNa"_IGEref_out01a_filtered_genetagged.tsv" > $SNa"_IGEref_out02_sepparated.tsv" && rm $SNa"_IGEref_out01a_filtered_genetagged.tsv" && rm $SNa"_IGEref_out01_genetagged.tsv"
	# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
	# and create column with BASE readname (no :A or :B)
		awk 'BEGIN {OFS = "\t"} {
		a = $5
		gsub(/_LTR/,"",a)
		c = substr($4, 1, length($4)-2)
		print $0"\t"a"\t"c
		}' < $SNa"_IGEref_out02_sepparated.tsv" > $SNa"_IGEref_out03pre_TEbase.tsv" && rm $SNa"_IGEref_out02_sepparated.tsv"
	wc_predup=$(wc -l $SNa"_IGEref_out03pre_TEbase.tsv" | awk '{print $1}')
	# remove duplicate reads (where both :A and :B version of the same reads were picked up
	awk '!seen[$21]++' $SNa"_IGEref_out03pre_TEbase.tsv" > $SNa"_IGEref_out03_TEbase.tsv" && rm $SNa"_IGEref_out03pre_TEbase.tsv"
	echo " ---- $(wc -l $SNa"_IGEref_out03_TEbase.tsv" | awk '{print $1}') reads ($wc_predup before removing duplicates)" >> $SNa"_PART3to5_"$logname".log"
	echo " <-- done processing PART1 output at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
}

combine_hits_of_each_TE()
{
	echo " --> start combining reads at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
	# create list of all TEs in dataset
	cut -f20 $1 | sort | uniq > $SNa"_IGEref_out03a_uniqueTEs.tsv"
	# to create collection file, first make sure this file does not yet exist
	rm -f $SNa"_IGEref_out12_OUTPUT.tsv"
	# loop through all TEs in the dataset
	while read TE
	do
		# grep TEs from combined data set
		awk -v t="$TE" '{if ($20 == t) print $0;}' $1 > "tmp."$SNa"_"$TE"_IGEref_out04.tsv"
		if [ -s "tmp."$SNa"_"$TE"_IGEref_out04.tsv" ]; then	
			###### SPLIT into TE-GENE and GENE-TE
			awk '{if ($7 == "TE-GENE") print $0;}' "tmp."$SNa"_"$TE"_IGEref_out04.tsv" > "tmp."$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv"
			# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
			# header: Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllBreakpoints(TE)
			if [ -s "tmp."$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_IGEref_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_IGEref_out06_TE-GENE.tsv"
			fi
			awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$SNa"_"$TE"_IGEref_out04.tsv" > "tmp."$SNa"_"$TE"_IGEref_out05_GENE-TE.tsv"
			if [ -s "tmp."$SNa"_"$TE"_IGEref_out06_GENE-TE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_IGEref_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,6,7,17,10 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse -d 20 > "tmp."$SNa"_"$TE"_IGEref_out06_GENE-TE.tsv"
			fi
			###### COMBINE TE-GENE and GENE-TE
			cat "tmp."$SNa"_"$TE"_IGEref_out06_"*".tsv" > "tmp."$SNa"_"$TE"_IGEref_out08.tsv"
			bedtools sort -i "tmp."$SNa"_"$TE"_IGEref_out08.tsv" > "tmp."$SNa"_"$TE"_IGEref_out09.tsv"
			# create list of all genes
			# first, replace "," with tab, then sort, then pick unique strings, then remove ""
			cut -f12 "tmp."$SNa"_"$TE"_IGEref_out09.tsv" | tr ',' '\n' | sort | uniq | sed '/\./d' > "tmp."$SNa"_"$TE"_IGEref_out11_genes.tsv"
			# this list of genes is used to add non-anchor hits to ANCHOR hits, if they overlap with the same gene and TE
			while read GENE
			do
				# create output file that contains for every gene all the insertions of that TE
				awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$13"\t"$7"\t"$8"\t""mRNA";}' "tmp."$SNa"_"$TE"_IGEref_out09.tsv" >> $SNa"_IGEref_out12_OUTPUT.tsv"
			done < "tmp."$SNa"_"$TE"_IGEref_out11_genes.tsv"
			rm -f "tmp."$SNa"_"$TE*
		fi
	done < $SNa"_IGEref_out03a_uniqueTEs.tsv"
	# remove empty lines
	grep . $SNa"_IGEref_out12_OUTPUT.tsv" > $SNa"_IGEref_out12a_OUTPUT.tsv"
	rm $SNa"_IGEref_out03a_uniqueTEs.tsv"
	rm $SNa"_IGEref_out03_TEbase.tsv"
	echo " <-- done combining reads at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
}

check_for_splicesites()
{
	echo " --> start checking for breaks at splice sites at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
	# now, every line is one incident of a gene-TE chimera
	while read line
	do
		# convert line to bed format
		# header: Chr(genome)\tStart(genome)[this is the breakpoint minus 1]\tEnd(genome)[this is the breakpoint]\tTE|Gene|TE-GENEorGENE-TE\t"."\tStrand(genome)
		echo $line | awk '{ print $2"\t"$4-1"\t"$4"\t"$5"|"$1"|"$6"\t.\t"$3}' > "tmp."$SNa"_IGEref_out13.bed"
		# $a is TE-GENE or GENE-TE	
		a=$(echo $line | awk '{print $6}')
		# $b is the gene
		b=$(echo $line | awk '{print $1}')
		# intersect with FEATURES file to get all ovelapping features (output is ;-separated)
		d=$(bedtools intersect -a "tmp."$SNa"_IGEref_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa -s | awk '{print $10}' | paste -sd ";" -)
		# check whether breakpoint on the genome overlaps with splice donor site (in the case of a GENE-TE fragment)
		# or with splice acceptor site (in the case of a TE-GENE fragment). the REF files _SPLICE_DONORS.bed and 
		# _SPLICE_ACCEPTORS.bed have a +/- 10nt, but the precise breakpoint is extracted here (which is the sequence
		# number of the last nucleotide that is still part of the exon for SPLICE_DONORS and the first exonic nt 
		# for SPLICE_ACCEPTORS.
		# for GENE-TE fragments, the exon is the splice donor
		if [ $a == "GENE-TE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_DONORS.bed" > "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$SNa"_IGEref_out13.bed" -b "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
		# for TE-GENE fragments, the exon is the splice acceptor
		elif [ $a == "TE-GENE" ]
		then
			grep $b $REFpath$REFbase"_SPLICE_ACCEPTORS.bed" > "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
			c=$(bedtools intersect -a "tmp."$SNa"_IGEref_out13.bed" -b "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed" -loj -wa | head -n 1 | awk '{print $10}' | cut -d "|" -f2)
			rm "tmp."$SNa"_IGEref_out14_ref_for_overlap.bed"
		else
			c="N/A"
		fi
		if [[ $c = "" ]]; then c="."; fi
		# the loop aboce created the unix variable $c. here, this variable is appended to the main file.
		echo $line | awk -v c="$c" -v d="$d" 'BEGIN {OFS=FS = "\t"} { print $0"\t"c"\t"d}'
		rm "tmp."$SNa"_IGEref_out13.bed"
	done < $1 > $SNa"_IGEref_out15.tsv"
	echo " <-- done checking for breaks at splice sites at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
}

split_TE_breakpoints()
{
	echo " --> start tyding up TE breakpoints at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
	# for better readability, the concatenated breakpoints on the TE are pooled here (with the number of occurrences in brackets)
	cat $1 | while read line; do echo $line | awk '{print $8}' | awk -v RS=',' '{print$0}' | sort | uniq -c | awk '{if (NR!=1) printf$2"("$1"),"}'| sed 's/,$//'; done > $SNa"_newcolb"
	# for final output, get rid of original TE-breakpoint field
	awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13;}' $SNa"_IGEref_out15.tsv" > $SNa"_newcola"
	# append pooled TE breakpoints to final output
	paste $SNa"_newcola" $SNa"_newcolb" > $SNa"_chimericreads_final.tsv"
	rm $SNa"_newcola" $SNa"_newcolb"
	echo " <-- done tyding up TE breakpoints at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
}


################################################################################
################################################################################
# execute:

# change to wd
cd $wd

# create .log file
logname=$(date | awk '{gsub(/\:/,"-",$5); print $4$3$2"_"$5}')
echo "====================" > $SNa"_PART3to5_"$logname".log"
echo "|| TEchim - PART3 || " >> $SNa"_PART3to5_"$logname".log"
echo "====================" >> $SNa"_PART3to5_"$logname".log"
echo "Parameters:" >> $SNa"_PART3to5_"$logname".log"
echo "Working directory:" "$wd/$SNa"_IGEref_CONTAINER"" >> $SNa"_"$logname".log"
echo "Location of PART1 output:" "$path_to_PART1_output" >> $SNa"_"$logname".log"
echo "Sample name:" "$SNa" >> $SNa"_PART3to5_"$logname".log"
echo "Reference files:" "$REFpath$REFbase" >> $SNa"_PART3to5_"$logname".log"
echo "--------------------------------" >> $SNa"_PART3to5_"$logname".log"
echo " --> starting at ... $(date)" >> $SNa"_PART3to5_"$logname".log"
echo "--------------------------------" >> $SNa"_PART3to5_"$logname".log"

find_matching_IGEs

create_IGE_reference $SNa"_IGEs.bed"

list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}')
for SNo in $list_of_snum
do
	list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2"\t"$3}' | grep $snum | awk '{print $2}')
	for LNo in $list_of_lanes
	do
		align_and_filter $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_out3_1.fasta" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_out3_2.fasta" &
	done
done
wait
	
list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}')
for SNo in $list_of_snum
do
	list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2"\t"$3}' | grep $snum | awk '{print $2}')
	for LNo in $list_of_lanes
	do
		blast_on_longreads $SNa"_IGEref_"$SNo"_"$LNo$"_out5_TExGENES.sam" $path_to_PART1_output$SNa"_"$SNo"_"$LNo"/"$SNa"_"$SNo"_"$LNo$"_LOOKUP.sorted.tsv.gz" &
	done
done
wait

list_of_snum=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2}' | awk '!seen[$0]++ {print $0}')
for SNo in $list_of_snum
do
	list_of_lanes=$(find $path_to_PART1_output -maxdepth 1 -name "$SNa"_S"*" | rev | cut -d "/" -f 1 | rev | awk '{gsub(/_/,"\t"); print $2"\t"$3}' | grep $snum | awk '{print $2}')
	for LNo in $list_of_lanes
	do
		create_summary_table $SNa"_IGEref_"$SNo"_"$LNo$"_out9_TExGENES_blastedreads.tsv" &
	done
done
wait

cd $wd
# create out10_combined.sorted.bed
cat $wd"/IGE_"$SNa*"/"$SNa*"_out10_breakpoints.bed" | bedtools sort -i - > $SNa"_IGEref_in10_combined.sorted.bed"

# add a gene tag to each read (i.e. genomic location). this is done using a gene-only version of the gtf file.
bedtools intersect -wa -a $SNa"_IGEref_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj -s > $SNa"_IGEref_out01_genetagged.tsv"

process_PART1_output $SNa"_IGEref_out01_genetagged.tsv"

combine_hits_of_each_TE $SNa"_IGEref_out03_TEbase.tsv"

check_for_splicesites $SNa"_IGEref_out12a_OUTPUT.tsv"

split_TE_breakpoints $SNa"_IGEref_out15.tsv"

echo " <-- all done at ... $(date)" >> $SNa"_PART3to5_"$logname".log"

