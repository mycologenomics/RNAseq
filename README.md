# RNAseq
Alignment of RNAseq reads to reference genome using hisat2

Hisat2 is a splice-aware aligner, and has been shown to be the gold standard for RNAseq alignment and has superseded Tophat 2. Hisat2 v2.0.4 is already on the HPC, but you can update to v2.2.1 (which I'm using) and install it via anaconda 

Firstly, you need to prepare your reference genome, so you need your reference genome in fasta format and a gtf. As an example, I am using the latest version of the Candida auris B8441 genome.

>source activate hisat2_env

Extract splice sites and exons from the gtf, which will be used when making the reference genome indexes for alignments - this is how hisat2 is splice aware:

>hisat2_extract_splice_sites.py GCA_002759435.2_Cand_auris_B8441_V2_genomic.gtf > GCA_002759435.2_Cand_auris_B8441_V2_genomic.ss
>hisat2_extract_exons.py GCA_002759435.2_Cand_auris_B8441_V2_genomic.gtf > GCA_002759435.2_Cand_auris_B8441_V2_genomic.exon

Now build the HFM index:

hisat2-build -p 16 --ss GCA_002759435.2_Cand_auris_B8441_V2_genomic.ss --exon GCA_002759435.2_Cand_auris_B8441_V2_genomic.exon GCA_002759435.2_Cand_auris_B8441_V2_genomic.fa genome

This will create multiple files with the prefix 'genome'

Now to align!

>hisat2 -p 8 --rg-id=VH00344_2_AAAKMY2M5 --rg SM:PR_0h_REP1 --rg LB:rnaseq --rg PL:ILLUMINA --rg PU:ATCACG -x genome --dta --rna-strandness RF -1 R1.fastq.gz -2 R2.fastq.gz -S PR_0h_REP1.sam

- ‘-p 8’ tells HISAT2 to use eight CPUs for bowtie alignments.
- ’–rna-strandness RF’ specifies strandness of RNAseq library. We will specify RF since the TruSeq strand-specific library was used to make these libraries. See here for options.
- ’–rg-id $ID’ specifies a read group ID that is a unique identifier.
- ’–rg SM:$SAMPLE_NAME’ specifies a read group sample name. This together with rg-id will allow you to determine which reads came from which sample in the merged bam later on.
- ’–rg LB:$LIBRARY_NAME’ specifies a read group library name. This together with rg-id will allow you to determine which reads came from which library in the merged bam later on.
- ’–rg PL:ILLUMINA’ specifies a read group sequencing platform.
- ’–rg PU:$PLATFORM_UNIT’ specifies a read group sequencing platform unit. Typically this consists of FLOWCELL-BARCODE.LANE
- ’–dta’ Reports alignments tailored for transcript assemblers.
- ‘-x /path/to/hisat2/index’ The HISAT2 index filename prefix (minus the trailing .X.ht2) built earlier including splice sites and exons.
- ‘-1 /path/to/read1.fastq.gz’ The read 1 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
- ‘-2 /path/to/read2.fastq.gz’ The read 2 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
- ‘-S /path/to/output.sam’ The output SAM format text file of alignments.

Now convert to BAM: 

>samtools view -u PR_0h_REP1.sam | samtools sort -@ 8 - PR_0h_REP1

This BAM file can be used as input for Stringtie to generate expression estimates in 'reference only' mode. This reduced run time. However, Stringtie can predict the transcripts present in each library by dropping the -G option and not providing the gtf/gff3.

stringtie --rf -p 28 -B -j 3 -c 7.5 -g 1 -v -G GCA_002759435.2_Cand_auris_B8441_V2_genomic.no_genes.gtf -o PR_0h_REP1.transcripts.gtf -A PR_0h_REP1.abundances.csv PR_0h_REP1.bam

Where: --rf tells StringTie that our data is stranded and to use the correct strand specific mode
-o output file name for assembled transcripts
-A output file name for gene abundance estimates
-p is the number of CPUs to use
-j
-c
-v verbose output
-G gtf/gff3 file
-g
-B create Ballgown specific output

The Ballgown specific output option creates ctab files, which are required for Ballgown. Make a csv. For example, I wanted to compare timepoints 0h and 6h, so my file (pheno_0h_vs_6h.csv) looked like this:

ids,type,path
0h_REP1,0h,extdata/0h_REP1
0h_REP2,0h,extdata/0h_REP2
0h_REP3,0h,extdata/0h_REP3
6h_REP1,6h,extdata/6h_REP1
6h_REP2,6h,extdata/6h_REP2
6h_REP3,6h,extdata/6h_REP3

The Ballgown analysis can be completed in R:

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Load phenotype data from a file we saved in the current working directory
pheno_data_PR = read.csv("pheno_0h_vs_6h.csv")

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data_PR$path), pData=pheno_data_PR)

# Display a description of this object
bg

# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Save the ballgown object to a file for later use
save(bg, file='bg_PR.rda')

# Perform differential expression (DE) analysis with no filtering
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts, "PR0h_vs_6h_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "PR0h_vs_6h_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_transcripts, "PR0h_vs_6h_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "PR0h_vs_6h_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

# Output the signifant gene results to a pair of tab delimited files
write.table(sig_transcripts, "PR0h_vs_6h_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes, "PR0h_vs_6h_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
