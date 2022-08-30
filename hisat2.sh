#!/bin/sh

## Script to align reads to reference using splice aware hisat2
## Input files required
## $1 --> 1st of paired read in fastq format
## $2 --> 2nd of paired read in fastq format
## $3 --> output prefix
## Versions used:
## hisat2 v2.2.1

mkdir -p -v /rds/general/project/fisher-candida-results/ephemeral/$3

ref_dir=/rds/general/project/fisher-candida-reference/live/
ref=$ref_dir/GCA_genome
output=/rds/general/project/fisher-candida-results/ephemeral/$3
results_dir=/rds/general/project/fisher-candida-results/live/RNAseq/$3
rgid=VH00344_2_AAAKMY2M5
rgpu=$4

cp -v $1 $output/R1.fastq.gz
cp -v $2 $output/R2.fastq.gz

hisat2 -p 8 --rg-id=$rgid --rg SM:$3 --rg LB:rnaseq --rg PL:ILLUMINA --rg PU:$rgpu -x $ref --dta --rna-strandness RF -1 $output/R1.fastq.gz -2 $output/R2.fastq.gz -S $output/$3.sam

samtools view -u $output/$3.sam | samtools sort -@ 8 - $output/$3

stringtie --rf -p 28 -j 3 -c 7.5 -g 1 -v -G $ref_dir/GCA_002759435.2_Cand_auris_B8441_V2_genomic.no_genes.gtf -e -o $output/$3.transcripts.gtf -A $output/$3.abundances.csv $output/$3.bam

cp -v $output/$3.bam $results_dir
cp -v $output/$3.transcripts.gtf $results_dir
cp -v $output/$3.abundances.csv $results_dir
