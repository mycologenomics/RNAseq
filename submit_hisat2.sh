#!/bin/sh
 
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb
## This tells the batch manager to re-run job with parameter varying from 1 to N in steps on step-size
#PBS -J 1-12

## OTHER OPTIONAL PBS DIRECTIVES

module load samtools
module load stringtie
module load anaconda3/personal
source activate hisat2_env

/rds/general/project/fisher-candida-analysis/live/hisat2.sh $(head -$PBS_ARRAY_INDEX /rds/general/project/fisher-candida-analysis/live/arglist_hisat2.txt | tail -1)
