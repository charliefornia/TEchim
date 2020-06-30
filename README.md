# TEchim

TEchim is an analysis pipeline with 6 key functions:

1. Generation of the required support files
2. Detection of transposon-gene breakpoint-spanning reads in gDNA data
3. Detection of transposon-gene breakpoint-spanning reads in cDNA data
4. Generation of immobile genetic element (IGE) analysis
5. Quantification of LTR-gene spanning reads
6. Quantification of locus-specific breakpoint-spanning reads

------------------------------------------------------------------------

TEchim – function 2/3
part 1 – on every sample/lane separately
1.) merge reads
If 2 fastq files exist, then these reads are merged using FLASH. The command for this step is:
flash original_1.fastq original_2.fastq -z [compresses the output file] -x 0.15 [maximum allowed ratio between the number of mismatched base pairs and the overlap length] -M 170 [maximum overlap] -o sample_out1 -t number_of_cores -q [quiet]

The merged reads are in the file sample_out1.extendedFrags.fastq.gz. The file sample_out1.notCombined_1.fastq.gz can be added to the merged reads, because it has the same orientation as the merged reads.

2.) create fasta
Next, in-silico paired-end reads are generated. A paired-end read can only be generated from a complete contiguous read. The resulting _1.fasta will represent the coding strand, regardless how the input was stranded.

Only input reads are selected that are larger than 2x the in-silico fasta read length, that is set in the parameters (fasta_length). Hence, this value should be chosen accordingly. Longer fasta_length values allow more precise mapping, but reduce the total number of reads that can be used to map chimeric transcripts.

-	First, the sample_out1.combined.fastq.gz and sample_out1.notCombined_1.fastq.gz files are concatenated, then split up into 4 separate files, using awk. (a) the first line contains the readname. In the reads from file original_1.fastq, an ":A" is added to the end of the readname (using sed), in order not to confuse them with (otherwise identical) readnames from original_2.fastq. (b) the line containing the sequence is cropped to fasta_length (taking the FIRST n nucleotides) (c) the line with a + is kept (d) the quality scores are also cropped to fasta_length
    -> these files are used to generate in_silico_fasta_1.fasta if the first read of the input file represented the coding strand, and _2 if input read 2 was coding. 
-	Second, the (b) sequence line is reversed and cropped from beginning (i.e. original end) and the reverse complement strand is generated. The (d) quality score is also reversed and cropped from the beginning (i.e. original end).
    -> these files are used to generate in_silico_fasta_2/1.fastq
-	Third, the procedure above is repeated for the sample_out1.notCombined_2.fastq.gz file, and the reads are appended to the in_silico_fasta_X.fastq files.
-	Fourth, a lookup table is generated, where each line contains the read name and a long, contiguous sequence (the one that was used to generate in_silico_fasta_X.fastq). If the input library was stranded, then the sequence is representing the mRNA sequence (e.g. ATG for start codon, and AAAAA(N) for the poly-A signal). This depends on the parameter stranded. If original_1.fastq was the coding strand, then stranded should be set to 1. If stranded=1, then the reads taken from sample_out1.combined.fastq.gz and sample_out1.notCombined_1.fastq.gz stay unchanged, and sample_out1.notCombined_2.fastq.gz is reverse-complemented and appended. If stranded=2, then the opposite is done. The file sample_LOOKUP.tsv is then sorted and compressed.

3.) align and filter
Next, three streams are performed to identify chimeric molecules.

STREAM 1:  the in-silico reads are first mapped to a masked reference genome. The STAR aligner is used, and detection of chimeric alignments is switched on. The comandline is:
STAR --runThreadN $nc --genomeDir $REFpath"STAR_"$REFbase --readFilesIn $1 $2 --chimSegmentMin 20 --chimOutType WithinBAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $SNa"_S"$SNo"_L"$LNo$"_out4_"
Read names of pairs that span gene-transposon breakpoints are extracted, and the matching long reads are blasted to (A) transposon-free reference genome and (B) a transposon-only “genome”. BLAST identifies the precise breakpoints. 

STREAM 2: STAR aligner also identifies chimeric transcripts. To extract those that span a transposon-gene breakpoint,  four samflag filters are applied:
-	samtools view -f 65 [1-read paired; 64-first in pair] -F 48 [NOT: 16-read reverse; 32-mate reverse]
-	samtools view -f 97 [1-read paired; 32-mate reverse; 64-first in pair] -F 16 [NOT: 16-read reverse]
-	samtools view -f 81 [1-read paired; 16-read reverse; 64-first in pair] -F 32 [NOT: 32-mate reverse]
-	samtools view -f 113 [1-read paired; 16-read reverse; 32-mate reverse; 64-first in pair]
Depending on the flags, and whether the first- or second read maps to a transposon, each read is either labelled TE-GENE or GENE-TE, the strand of both the genomic- and transposon- site are determined, and the breakpoints in gene and transposons are extracted. Since these “hits” are not based on BLAST of contiguous reads, the breakpoints are subject to a range. For the range, the maximum fragment size of the FLASH merge step is taken. The hits from STREAM 2 are only used if the long-read version has a positive hit in STREAM 1 for the genomic section (but not TE was found using BLAST).

STREAM 3: STAR is also performed on the long-reads directly. This alignment also gives precise breakpoints, but work at lot less efficiently than the other 2 streams. Nevertheless, reads that have NOT been identified as evidence for chimera in streams 1 and 2 are added to the final output file.

The final output file contains a line for every sequencing reads that spans a transposon-gene breakpoint. The format is:

2L	123500	123546	SRR070434.3715387.1:A|QUASIMODO2_LTR|plus|TE-GENE|S1|L1|144-294|0	.	-

Part1 is run on every lane and sample.

Part 2 – combine samples/lanes
1.) process_P1out_TE
First, all existing output files from part1 in the destination folder are combined (concatenated). The bedfile is intersected with all annotated “GENES” in the reference genome, and a gene tag is added to each read (or -1, if the genomic section of the read does not overlap with an annotated gene).

Next, duplicate reads are removed (i.e. reads where both the :A and :B version were called). 

2.) combine_hits_of_each_TE
A list of all TEs is first created. For TE_LTR, only the basename is used. Then, the data is looped for each transposon. In the loop, these steps are run:
-	The output file is split into TE-GENE and GENE-TE reads.
-	bedtools merge is run on both these files, merging reads that map within 20 nucleotides of each other. 
The output format is

Chr(genome)|Start(genome)|End(genome)|TE|"."|Strand(genome)|NumberOfSamples|NumberOfReads|Breakpoint(genome)|Strand(TE)|TE-GENEorGENE-TE|GeneNames|AllSamples|AllBreakpoints(TE)

Chr (genome)	Taken from bedtools merge
Start (genome)	Taken from bedtools merge
End (genome)	Taken from bedtools merge
TE	[distinct], TE name. If LTR is also detected it is printed here 
.	[mode] Required for bedfile
Strand (genome)	[mode] Most frequent genome strand is taken
Number of samples	[count_distinct] Number of different S’ (e.g. S1, S1, S2, S5 -> 3)
Number of reads	[count_distinct] Unique number of reads (that don’t have :A/:B)
Breakpoint (genome)	[mode] Most frequent breakpoint is taken
Orientation (TE)	[mode] Most frequent orientation is taken
TE-GENE or GENE-TE	[mode] This will be the same for all reads in sub-section
Gene names	[distinct] All gene names are taken
All samples	[collapse] List of all samples that contributed reads
All breakpoints	[collapse] List of all TE breakpoints or range-of-breakage

-	Both files are then combined again and sorted. 
-	Next, all genes are taken and looped through for each transposon. Because of this step, some distinct transposon-gene reads might sometimes be counted twice. This will happen, when the breakpoint is located inside a gene and an intron of a larger, surrounding gene.
In the gene loop, the output file is re-structured

Gene|Chr|Strand(genome)|Breakpoint(genome)|TE|TE-GENE or GENE-TE|orientation (te)|Number of samples|Number of reads|assay|AllSamples|AllBreakpoints

After all TEs are looped through, the empty lines are removed.

3.) check_for_TE_splicesites

Next, every line in output file is checked for proximity to a splice site of the according gene. First, the line is converted to bedfile format. Then, this bedfile is intersected with $REFpath$REFbase”_FEATURES.bed”, which contains all features of every gene (e.g. intron, UTR, exon, etc.). The entire list of overlapping features is appended to each line.
After this, TEchim checks whether the breakpoint overlaps with a splice junction of the host gene. For TE-GENE fragments, splice acceptor sites are overlapped, and for GENE-TE fragments, splice donor sites are taken.
Finally, all fields are pasted into a tab-separated output file.

4.) split_TE_breakpoints
This step splits and pools the TE breakpoints for each line.

5.) add_expression_levels
The last segment quantifies the number of transcripts with- and without transposon insertions. This analysis depends on whether TE-GENE, or GENE-TE fragments were detected. First, 


two regions are determined for every line in the output file.
For TE-GENE fragments, region A stretches from the start of the gene up until the breakpoint and region B starts at the breakpoint and extends downstream of the gene for a length equivalent to the maximum fragment length (as determined after the fastq-merging step in PART 1).
For GENE-TE fragments, region A stretches from the breakpoint up until the end of the gene and region B ends at the breakpoint and extends upstream of the gene for the maximum fragment length.
Two values are then calculated for every sequencing lane and sample. One is the number of properly paired reads where one mate maps onto region B, and the other mate onto the inserted transposon. And the other value is the number of properly paired reads where one mate maps onto region B, and the other mate onto region A.
All these steps are strand-sensitive, and only count reads that match the according strand of the gene.


TEchim function 4 - IGE control analysis

1.) split_CDS

The program starts in the same working directory as Parts 1 & 2, and creates a folder called “IGE_COLLECTION”. Next, the top 10% of coding sequences (CDS) are counted, and the CDS are divided into 10 randomly composed groups.

2.) calculate_TE_coverage

Next, the expression levels of both the positive- and the negative strand of all tested transposons in all samples is calculated. This is done in a strand-specific manner.

3.) find_matching_IGEs

Based on the average expression levels of all transposons, 10 sets of IGEs with matching expression levels are picked.

4.) create_IGE_reference

10 new reference genomes are created, where each set of IGEs is first removed, and then added as separate chromosomes (named IGE_...), analogous to transposons in the original reference genome adaptations.

The next steps are performed identical to the analysis of Transposons, for each of the 10 IGE sets.

