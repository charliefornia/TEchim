#!/bin/bash

################################################################################
# TITLE: TEchim - part 2
# VERSION: 0.1.1 (dev)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 09/05/2019 (dd/mm/yyyy)
# DESCRIPTION: tbd
################################################################################

################################################################################
# REQUIREMENTS:
# - bedtools
################################################################################

# BUILD GENE FILE
# awk 'BEGIN {OFS = "\t"} {a = $0; b = $10; gsub(/"/,"",b); gsub(/;/,"",b) ; if ($3 == "gene") {print $1"\t"$4"\t"$5"\t"b"\t"$6"\t"$7}}' dmel625_v03.gtf | bedtools sort -i - > dmel625_v04_ONLYgenesFROMgtf.bed

# set parameters
wd=/Users/CT/Dropbox/CloudDesktop/TEchim_cloud/testing
input=/Users/CT/Dropbox/CloudDesktop/TEchim_cloud/2018MARCH_TEchim_combined.sorted.bed
gene_ref=/Users/CT/Dropbox/CloudDesktop/REF_cloud/dmel625_v04_ONLYgenesFROMgtf.bed
samplename=2018MARCH_TEchim

# filters for "anchor" insertion
min_repl=2
min_reads=4

# change to wd
cd $wd

# add a gene tag to each read (i.e. genomic location). this is done using a gene-only version of the gtf file.
bedtools intersect -wa -a $input -b $gene_ref -loj -s > $samplename"_out01_genetagged.tsv"

# separate | delimited field
sed -e $'s/|/\t/g' $samplename"_out01_genetagged.tsv" > $samplename"_out02_sepparated.tsv"

# generate column that contains the "basic" TE name i.e. TE_LTR ==> TE
awk 'BEGIN {OFS = "\t"} {
	a = $5
	gsub(/_LTR/,"",a)
	c = substr($4, 1, length($4)-2)
	print $0"\t"a"\t"c
	}' < $samplename"_out02_sepparated.tsv" > $samplename"_out03_TEbase.tsv"

# create list of all TEs in dataset
cut -f20 $samplename"_out03_TEbase.tsv" | sort | uniq > $samplename"_out03a_uniqueTEs.tsv"

# to create collection file, first make sure this file does not yet exist
rm -f $samplename"_out06.maininsertions.tsv"
rm -f $samplename"_out12_OUTPUT.tsv"

# loop through all TEs in the dataset
while read TE
do
	# grep TEs from combined data set
	awk -v t="$TE" '{if ($20 == t) print $0;}' $samplename"_out03_TEbase.tsv" > "tmp."$samplename"_"$TE"_out04.tsv"
	
	###### SPLIT into TE-GENE and GENE-TE ?!
	awk '{if ($7 == "TE-GENE") print $0;}' "tmp."$samplename"_"$TE"_out04.tsv" > "tmp."$samplename"_"$TE"_out05_TE-GENE.tsv"
	awk '{if ($7 == "GENE-TE") print $0;}' "tmp."$samplename"_"$TE"_out04.tsv" > "tmp."$samplename"_"$TE"_out05_GENE-TE.tsv"
	
	# merge, counting unique read NAMES, and allowing for 20nt range of precise insertion site
	# header: 
	bedtools merge -i "tmp."$samplename"_"$TE"_out05_TE-GENE.tsv" -c 5,12,13,8,21,3,10,6,7,17 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,mode,distinct -d 20 > "tmp."$samplename"_"$TE"_out06_TE-GENE.tsv"
	bedtools merge -i "tmp."$samplename"_"$TE"_out05_GENE-TE.tsv" -c 5,12,13,8,21,3,10,6,7,17 -o mode,mode,mode,count_distinct,count_distinct,mode,mode,mode,mode,distinct -d 20 > "tmp."$samplename"_"$TE"_out06_GENE-TE.tsv"

########### FOR ANCHOR
	# filter to keep only insertions that were detected in at least $min_repl separate biological replicates with at least $min_reads reads
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {print a}}' < "tmp."$samplename"_"$TE"_out06_TE-GENE.tsv" > "tmp."$samplename"_"$TE"_out07_TE-GENE.tsv"
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {print a}}' < "tmp."$samplename"_"$TE"_out06_GENE-TE.tsv" > "tmp."$samplename"_"$TE"_out07_GENE-TE.tsv"
	# combine TE-GENE and GENE-TE, then sort
	cat "tmp."$samplename"_"$TE"_out07_TE-GENE.tsv" "tmp."$samplename"_"$TE"_out07_GENE-TE.tsv" > "tmp."$samplename"_"$TE"_out08.tsv"
	bedtools sort -i "tmp."$samplename"_"$TE"_out08.tsv" > "tmp."$samplename"_"$TE"_out09.tsv"

########### FOR rest
	# filter to keep only insertions that were detected in at least $min_repl separate biological replicates with at least $min_reads reads
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {} else {print a}}' < "tmp."$samplename"_"$TE"_out06_TE-GENE.tsv" > "tmp."$samplename"_"$TE"_out07b_TE-GENE_REST.tsv"
	awk -v mreps="$min_repl" -v mreads="$min_reads" 'BEGIN {OFS = "\t"} {a = $0 ; if ( $7 >= mreps && $8 >= mreads) {} else {print a}}' < "tmp."$samplename"_"$TE"_out06_GENE-TE.tsv" > "tmp."$samplename"_"$TE"_out07b_GENE-TE_REST.tsv"
	# combine TE-GENE and GENE-TE, then sort
	cat "tmp."$samplename"_"$TE"_out07b_TE-GENE_REST.tsv" "tmp."$samplename"_"$TE"_out07b_GENE-TE_REST.tsv" > "tmp."$samplename"_"$TE"_out08b_REST.tsv"
	bedtools sort -i "tmp."$samplename"_"$TE"_out08b_REST.tsv" > "tmp."$samplename"_"$TE"_out09b_REST.tsv"


	# if there are several genes then split them up to separate columns
	# first, replace "," with tab
	sed -e $'s/,/\t/g' "tmp."$samplename"_"$TE"_out09.tsv" > "tmp."$samplename"_"$TE"_out10.tsv"
	cut -f13- "tmp."$samplename"_"$TE"_out10.tsv" | tr '\t' '\n' | sort | uniq | sed '/\./d' > "tmp."$samplename"_"$TE"_out11_genes.tsv"
	
	while read GENE
		do
			# create output file that contains for every gene with at least one anchor insertion, all the insertions of that TE
			awk -v g="$GENE" '{if ($13 ~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$12"\t"$7"\t"$8"\t""YES";}' "tmp."$samplename"_"$TE"_out09.tsv" >> $samplename"_out12_OUTPUT.tsv"
			awk -v g="$GENE" '{if ($13 ~ g) print g"\t"$1"\t"$6"\t"$9"\t"$4"\t"$11"\t"$10"\t"$12"\t"$7"\t"$8"\t""NO";}' "tmp."$samplename"_"$TE"_out09b_REST.tsv" >> $samplename"_out12_OUTPUT.tsv"
		done < "tmp."$samplename"_"$TE"_out11_genes.tsv"
	
	rm "tmp."$samplename*
done < $samplename"_out03a_uniqueTEs.tsv"






