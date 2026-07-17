#!/bin/bash

# Resources:
#SBATCH --time=0-08:00:00  # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=5G
#SBATCH --partition=short

# Environment:
#SBATCH --export=NONE


################################################################################
# TITLE: TEchim - PART 2 (gDNA)
# VERSION: 0.1.0 (dev)
# AUTHOR: Christoph Treiber, University of Oxford
# DATE: 17/07/2026 (dd/mm/yyyy)
# DESCRIPTION: This script combines output from PART1 for genomic DNA (gDNA)
# samples. It is a variant of TEChim_PART2_forSLURM adapted for gDNA input:
# gDNA reads have no mRNA strand and do not undergo splicing, so this script
# drops the splice-donor/acceptor and splice-ratio analysis that only makes
# sense for cDNA/mRNA data (see TEChim_PART2_forSLURM for that pipeline).
# Genomic breakpoints are still annotated against gene features (exon, intron,
# UTR, CDS, etc.) from _FEATURES.bed.
################################################################################

################################################################################
# REQUIREMENTS:
# - STAR (https://github.com/alexdobin/STAR)
# - samtools
# - bedtools
# - blast
################################################################################
################################################################################
# load modules

module load rna-star
module load blast
module load parallel

################################################################################
################################################################################
# set parameters
wd=/path/to/wd					# working directory (no trailing "/")
path_to_PART1_output=$wd"/"		# with trailing "/"
SNa="EXP_NAME"					# Experiment name
nc=${SLURM_CPUS_PER_TASK:-1}	# number of cores (default:1)
REFpath=				# if left empty (recommended) then REFpath from PART1 is used
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
	echo "=====================" > $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "|| TEchim - PART 2 (gDNA) || " >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "=====================" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "Parameters:" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "Working directory:" "$wd" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "Location of PART1 output:" "$path_to_PART1_output" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "Sample name:" "$SNa" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "Reference files:" "$REFpath$REFbase" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "--------------------------------" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo " --> starting at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo "--------------------------------" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
}

process_P1out_TE()
{
	cd $wd
	echo " --> start processing PART1 output for TE at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	cat $path_to_PART1_output$SNa*"/"$SNa*"_out13_breakpoints.bed" | bedtools sort -i - > $SNa"_in10_combined.sorted.bed"
	# gDNA has no mRNA strand, so genes are tagged without requiring
	# strand agreement between the read and the annotated gene
	bedtools intersect -wa -a $SNa"_in10_combined.sorted.bed" -b $REFpath$REFbase"_GENES.bed" -loj > $SNa"_out01_genetagged.tsv"

	# separate | delimited field
	tr '|' '\t' < $SNa"_out01_genetagged.tsv" > $SNa"_out02_sepparated.tsv" && rm $SNa"_out01_genetagged.tsv"
	# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
	# and create column with BASE readname (no :A or :B)
		awk 'BEGIN {OFS = "\t"} {
		a = $5
		gsub(/_LTR/,"",a)
		c = substr($4, 1, length($4)-2)
		print $0"\t"a"\t"c
		}' < $SNa"_out02_sepparated.tsv" > $SNa"_out03pre_TEbase.tsv" && rm $SNa"_out02_sepparated.tsv"
	wc_predup=$(wc -l $SNa"_out03pre_TEbase.tsv" | awk '{print $1}')
	# remove duplicate reads (where both :A and :B version of the same reads were picked up
	awk '!seen[$21]++' $SNa"_out03pre_TEbase.tsv" > $SNa"_out03_TEbase.tsv" && rm $SNa"_out03pre_TEbase.tsv"
	rm $SNa"_in10_combined.sorted.bed"
	echo " ---- $(wc -l $SNa"_out03_TEbase.tsv" | awk '{print $1}') TE reads ($wc_predup before removing duplicates)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	echo " <-- done processing PART1 output for TE at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
}

combine_hits_of_each_TE()
{
	cd $wd
	echo " --> start combining TE reads at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	# create list of all TEs in dataset
	cut -f20 $1 | sort | uniq > $SNa"_out03a_uniqueTEs.tsv"
	# to create collection file, first make sure this file does not yet exist
	rm -f $SNa"_out12_OUTPUT.tsv"
	# loop through all TEs in the dataset
	while read TE
	do
		# grep TEs from combined data set
		awk -v t="$TE" '{if ($20 == t) print $0;}' $1 > "tmp."$SNa"_"$TE"_out04.tsv"
		if [ -s "tmp."$SNa"_"$TE"_out04.tsv" ]; then
			###### SPLIT into TE-GENE and GENE-TE
			awk '{if ($7 == "TE-GENE") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv"
			# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
			# header: Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllSamples|AllBreakpoints(TE)
			if [ -s "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,6,7,17,8,10,5 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse,collapse,distinct -d 20 > "tmp."$SNa"_"$TE"_out06_TE-GENE.tsv"
			fi
			awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$SNa"_"$TE"_out04.tsv" > "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv"
			if [ -s "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv" ]; then
				bedtools merge -i "tmp."$SNa"_"$TE"_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,6,7,17,8,10,5 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,distinct,collapse,collapse,distinct -d 20 > "tmp."$SNa"_"$TE"_out06_GENE-TE.tsv"
			fi
			###### COMBINE TE-GENE and GENE-TE
			cat "tmp."$SNa"_"$TE"_out06_"*".tsv" > "tmp."$SNa"_"$TE"_out08.tsv"
			bedtools sort -i "tmp."$SNa"_"$TE"_out08.tsv" > "tmp."$SNa"_"$TE"_out09.tsv"
			# create list of all genes
			# first, replace "," with tab, then sort, then pick unique strings, then remove ""
			cut -f12 "tmp."$SNa"_"$TE"_out09.tsv" | tr ',' '\n' | sort | uniq | sed '/\./d' > "tmp."$SNa"_"$TE"_out11_genes.tsv"
			# this list of genes is used to add non-anchor hits to ANCHOR hits, if they overlap with the same gene and TE
			while read GENE
			do
				# create output file that contains for every gene all the insertions of that TE
				# "assay" column is labelled gDNA so downstream tables are self-documenting
				awk -v g="$GENE" '{if ($12 ~ g && $4 !~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$7"\t"$8"\tgDNA\t"$13"\t"$14"\t"$15;}' "tmp."$SNa"_"$TE"_out09.tsv" >> $SNa"_out12_OUTPUT.tsv"
			done < "tmp."$SNa"_"$TE"_out11_genes.tsv"
			rm -f "tmp."$SNa"_"$TE*
		fi
	done < $SNa"_out03a_uniqueTEs.tsv"
	# remove empty lines
	grep . $SNa"_out12_OUTPUT.tsv" > $SNa"_out12a_OUTPUT.tsv"
	rm $SNa"_out12_OUTPUT.tsv"
	rm $SNa"_out03a_uniqueTEs.tsv"
	rm $SNa"_out03_TEbase.tsv"
	echo " <-- done combining TE reads at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
}

check_for_TE_features()
{
	cd $wd
	echo " --> start annotating TE breakpoints with gene features at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
	# for gDNA there is no splicing, so (unlike TEChim_PART2_forSLURM) this step
	# does not check for proximity to splice donor/acceptor sites. It only
	# annotates which gene feature (exon, intron, UTR, CDS, etc.) the genomic
	# breakpoint overlaps, using $REFpath$REFbase"_FEATURES.bed"
	while IFS=  read -r line
	do
		# convert line to bed format
		# header: Chr(genome)\tStart(genome)[this is the breakpoint minus 1]\tEnd(genome)[this is the breakpoint]\tTE|Gene|TE-GENEorGENE-TE\t"."\tStrand(genome)
		awk 'BEGIN{OFS="\t"}{print $2, $4-1, $4, $5"|"$1"|"$6, ".", $3}' <<< "$line" > "tmp."$SNa"_out13.bed"
		# intersect with FEATURES file to get all overlapping features (output is ;-separated)
		# gDNA has no mRNA strand, so this is not run with -s (strand-matching)
		d=$(bedtools intersect -a "tmp."$SNa"_out13.bed" -b $REFpath$REFbase"_FEATURES.bed" -loj -wa | awk '{print $10}' | paste -sd ";" -)
		# "N/A" placeholder keeps the column layout identical to the cDNA/mRNA
		# pipeline's splice-site column, so downstream postproc scripts
		# (e.g. filterPart2Output) work unmodified on either output
		c="N/A"
		awk -v c="$c" -v d="$d" 'BEGIN{OFS=FS="\t"}{print $0, c, d}' <<< "$line"
		rm "tmp."$SNa"_out13.bed"
	done < $1 > $SNa"_out15.tsv"
	rm $1
	echo " <-- done annotating TE breakpoints with gene features at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
}

split_TE_breakpoints()
{
        cd "$wd"
        echo " --> start tyding up TE breakpoints at ... $(date)" >> "$wd/$SNa_PART2gDNA_$logname.log"

        # pool TE breakpoints from col 12
        while IFS= read -r line; do
                awk -F'\t' -v OFS="\t" '{print $12}' <<< "$line" \
                | awk -v RS=',' '{print}' \
                | LC_ALL=C sort \
                | uniq -c \
                | awk '{if (NR!=1) printf "%s(%s),",$2,$1} END{print ""}' \
                | awk 'NF'
        done < "$1" > "${SNa}_newcolb"

        # pool TE breakpoints from col 11
        while IFS= read -r line; do
                awk -F'\t' -v OFS="\t" '{print $11}' <<< "$line" \
                | awk -v RS=',' '{print}' \
                | LC_ALL=C sort \
                | uniq -c \
                | awk '{if (NR!=1) printf "%s(%s),",$2,$1} END{print ""}' \
                | awk 'NF'
        done < "$1" > "${SNa}_newcolc"

        # drop original TE-breakpoint field(s) and keep the same column reordering
        awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$13,$14,$15}' "$1" > "${SNa}_newcola"

        # append pooled TE breakpoints
        paste "${SNa}_newcola" "${SNa}_newcolb" "${SNa}_newcolc" > "${SNa}_TE_chimericreads_final.tsv"

        rm -f "${SNa}_newcola" "${SNa}_newcolb" "${SNa}_newcolc"
        rm -f "$1"
        echo " <-- done tyding up TE breakpoints at ... $(date)" >> "$wd/$SNa_PART2gDNA_$logname.log"
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
check_for_TE_features $wd"/"$SNa"_out12a_OUTPUT.tsv"
split_TE_breakpoints $wd"/"$SNa"_out15.tsv"

echo " <-- all done at ... $(date)" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
echo "================================" >> $wd"/"$SNa"_PART2gDNA_"$logname".log"
