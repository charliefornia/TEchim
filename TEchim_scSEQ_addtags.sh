#!/bin/bash

################################################################################
# TITLE: Add cell barcode and UMI barcode to fasta files
# VERSION: 0.1.0 (beta)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 19/02/2020 (dd/mm/yyyy)
################################################################################

################################################################################
# This pipeline uses DropSeq tools
# (see https://github.com/broadinstitute/Drop-seq)
# Please cite:
# "Highly Parallel Genome-wide Expression Profiling of Individual Cells Using
# Nanoliter Droplets"; Evan Z. Macosko, ..., Steven A. McCarroll; Cell 2015
################################################################################

################################################################################
# REQUIREMENTS in $PATH:
# - DropSeqTools
# - java
################################################################################

usage () {
	echo "#########"
	echo "# Usage: TEchim_scSEQ_addtags.sh [options]"
	echo "#    -fq1 </PATH/TO/FILE/reads_with_barcodes.fastq(.gz)>"
	echo "#    -fq2 </PATH/TO/FILE/reads_mapping_to_genes.fastq(.gz)>"
	echo "#    -kit (1 - DropSeq, 2 - 10x-v2, 3 - 10x-v3 [default])"
	echo "#########"
}

while [ -n "$1" ]; do # while loop starts
    case "$1" in
    -fq1) 
    	FASTQ1="$2"
    	shift
    	;;
    -fq2)
        FASTQ2="$2"
        shift
        ;;
    -kit) 
    	kit="$2"
    	shift
    	;;
    *) echo "Option $1 not recognized" ;;
    esac
    shift
done

if [[ ! -f "$FASTQ1" ]] ; then
	echo "No FASTQ1 entered - use: -fq1 /PATH/TO/FILE"
	usage
	exit 1
fi
if [[ ! 
-f "$FASTQ2" ]] ; then
	echo "No FASTQ2 entered - use: -fq2 /PATH/TO/FILE"
	usage
	exit 1
fi
if [[ -z "$kit" ]] ; then
	kit="3"
fi

case "$kit" in
	"1")
		echo "using DropSeq settings"
		celltaglength=12
		UMIlength=8
		;;
	"2")
		celltaglength=16
		UMIlength=10
		;;
	"3")
		celltaglength=16
		UMIlength=12
		;;
	*) echo "KIT ID not valid"
		usage
		exit 1
		;;
esac

rand=$RANDOM

# get file basename
fbname=$(basename "$FASTQ2" | cut -d. -f1)

# get location of DropSeq tools. Dropseq tools must be in the $PATH.
pathtoDS=$(which Drop-seq_alignment.sh | rev | cut -d "/" -f2- | rev)

# STEP1: convert Fastq to SAM
java -jar $pathtoDS"/3rdParty/picard/picard.jar" FastqToSam F1=$FASTQ1 F2=$FASTQ2 O="./tmp."$rand"_01" SM=$fbname

# STEP2: add flag containing cell barcode
$pathtoDS"/TagBamWithReadSequenceExtended" SUMMARY="./tmp."$rand"_01to02_summary" BASE_RANGE="1-"$celltaglength BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 INPUT="./tmp."$rand"_01" OUTPUT="./tmp."$rand"_02"

# STEP3: add flag containing UMI barcode
$pathtoDS"/TagBamWithReadSequenceExtended" SUMMARY="./tmp."$rand"_02to03_summary" BASE_RANGE=$((celltaglength + 1))-$((UMIlength + celltaglength)) BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT="./tmp."$rand"_02" OUTPUT="./tmp."$rand"_03"

# STEP4: create new fastq.gz file where read name also contains cell barcode and UMI - separated with €
cat <(samtools view -H "./tmp."$rand"_03") <(samtools view "./tmp."$rand"_03" | awk '{OFS="\t"} {for (i=1;i<=NF;i++){if ($i ~/^XC:/) {xc=$i} else if ($i ~/^XM:/) {xm=$i}} {print $1"€"substr(xc,6)"€"substr(xm,6),$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}}') | samtools view -bh - | bedtools bamtofastq -i - -fq $fbname"_Cell_UMI_tagged.fastq"
gzip -f $fbname"_Cell_UMI_tagged.fastq"

rm "./tmp."$rand"_01"
rm "./tmp."$rand"_01to02_summary"
rm "./tmp."$rand"_02"
rm "./tmp."$rand"_02to03_summary"
rm "./tmp."$rand"_03"
