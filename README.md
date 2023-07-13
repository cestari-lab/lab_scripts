# lab_scripts
Set of scripts for Oxford nanopore sequencing data Alignment, RNA-seq and ChIP-seq:
Alignment and filtering - minimap and Samtools: Use this script for .fastq data alignment with minimap2 and obtain alignment statistics with Samtools. It uses samtools to convert .sam to .bam; sort and filter the alignments as required.
RNA-seq - featureCounts and edger: Use this script for counting mapped reads using package subread. Then analyze the counted data using edgeR to compare groups (e.g., control vs treatment) and obtain data statistics.
ChIP-seq - deeptools and macs3: Use this script for comparing chip-seq data (_sorted.bam files, e.g., chip vs input) to calculate enrichment of protein binding sites. It uses macs3 for peak calling and statistical analysis.

You can modify these scripts as you find necessary.
