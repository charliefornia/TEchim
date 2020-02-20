#!/bin/bash

################################################################################
# TITLE: Add cell barcode and UMI barcode to fasta files
# VERSION: 0.1.0 (beta)
# AUTHOR: Christoph Treiber, Waddell lab, University of Oxford
# DATE: 19/02/2020 (dd/mm/yyyy)
################################################################################

################################################################################
# This pipeline is largely based on the Macosko DropSeq Cookbook
# (see https://github.com/broadinstitute/Drop-seq)
# Please cite:
# "Highly Parallel Genome-wide Expression Profiling of Individual Cells Using
# Nanoliter Droplets"; Evan Z. Macosko, ..., Steven A. McCarroll; Cell 2015
################################################################################

################################################################################
# REQUIREMENTS:
# - DropSeqTools
# - java
################################################################################

################################################################################
# || VERY IMPORTANT :
# || Make sure that the cell- and UMI barcode lengths are correct !!!!
celltaglength=NUMBER_OF_READS
UMIlength=NUMBER_OF_READS
# || For 10x kit 3' - version2: celltaglength=16, UMIlength=10
# || For 10x kit 3' - version3: celltaglength=16, UMIlength=12
# || For Dropseq (Macosko, Chemgene beads): celltaglength=12, UMIlength=8
################################################################################

################################################################################
################################################################################
# set parameters
wd=$(pwd)						# working directory
SNa=NAME_OF_EXP					# sample name
FASTQ1=READS_1.fastq.gz			# input FASTQ file 1 (FULL PATH)
FASTQ2=READS_2.fastq.gz			# input FASTQ file 2 (FULL PATH)
SNo=1							# sample number (use integer from 1-n)
LNo=1							# sequencing lane number (use integer from 1-n)
################################################################################
################################################################################

cd $wd
# find location of DropSeq tools. Dropseq tools must be in the $PATH.
pathtoDS=$(which Drop-seq_alignment.sh | rev | cut -d "/" -f2- | rev)

# STEP1: convert Fastq to SAM
java -jar $pathtoDS"/3rdParty/picard/picard.jar" FastqToSam F1=$FASTQ1 F2=$FASTQ2 O="./tmp.01_"$SNa"_S"$SNo"_L"$LNo"_unmapped.bam" SM=$SNa

# STEP2: add flag containing cell barcode
$pathtoDS"/TagBamWithReadSequenceExtended" SUMMARY="./tmp.01to02_summary.txt" BASE_RANGE="1-"$celltaglength BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 INPUT="./tmp.01_"$SNa"_S"$SNo"_L"$LNo"_unmapped.bam" OUTPUT="./tmp.02_"$SNa"_S"$SNo"_L"$LNo"_unaligned_CellTagged.bam"

# STEP3: add flag containing UMI barcode
$pathtoDS"/TagBamWithReadSequenceExtended" SUMMARY="./tmp.02to03_summary.txt" BASE_RANGE=$((celltaglength + 1))-$((UMIlength + celltaglength)) BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT="./tmp.02_"$SNa"_S"$SNo"_L"$LNo"_unaligned_CellTagged.bam" OUTPUT="./tmp.03_"$SNa"_S"$SNo"_L"$LNo"_unaligned_CellUMITagged.bam"

# STEP4: create new fastq.gz file where read name also contains cell barcode and UMI - separated with €
cat <(samtools view -H "./tmp.03_"$SNa"_S"$SNo"_L"$LNo"_unaligned_CellUMITagged.bam") <(samtools view "./tmp.03_"$SNa"_S"$SNo"_L"$LNo"_unaligned_CellUMITagged.bam" | awk '{OFS="\t"} {for (i=1;i<=NF;i++){if ($i ~/^XC:/) {xc=$i} else if ($i ~/^XM:/) {xm=$i}} {print $1"€"substr(xc,6)"€"substr(xm,6),$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}}') | samtools view -bh - | bedtools bamtofastq -i - -fq "./"$SNa"_S"$SNo"_S"$LNo"_Cell_UMI_tagged.fastq"
gzip "./"$SNa"_S"$SNo"_S"$LNo"_Cell_UMI_tagged.fastq"

rm "./tmp.01_"$SNa"_S"$SNo"_L"$LNo"_unmapped.bam"
rm "./tmp.01to02_summary.txt"
rm "./tmp.02_"$SNa"_S"$SNo"_L"$LNo"_unaligned_CellTagged.bam"
rm "./tmp.02to03_summary.txt"
rm "./tmp.03_"$SNa"_S"$SNo"_L"$LNo"_unaligned_CellUMITagged.bam"


