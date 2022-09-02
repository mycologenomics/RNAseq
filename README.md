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

stringtie --rf -p 28 -j 3 -c 7.5 -g 1 -v -G GCA_002759435.2_Cand_auris_B8441_V2_genomic.no_genes.gtf -o PR_0h_REP1.transcripts.gtf -A PR_0h_REP1.abundances.csv PR_0h_REP1.bam

Where: --rf tells StringTie that our data is stranded and to use the correct strand specific mode
-o output file name for assembled transcripts
-A output file name for gene abundance estimates
-p is the number of CPUs to use
-j
-c
-v verbose output
-G gtf/gff3 file
-g
