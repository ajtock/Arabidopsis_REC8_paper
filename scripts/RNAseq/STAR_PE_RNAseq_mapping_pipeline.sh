#!/bin/bash

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Align paired-end RNA-seq reads to the
#  Arabidopsis TAIR10 reference genome assembly
# Software: STAR version 2.5.3a

# Example usage via condor submission system on hydrogen node7
# csmit -m 14G -c 24 "bash STAR_PE_RNAseq_mapping_pipeline.sh wt_RNAseq_Rep1"

i=$1

# align reads to reference genome using STAR
STAR --runThreadN 24 \
     --genomeDir /projects/ajt200/TAIR10/STAR_genome_index/ \
     --readFilesIn ./${i}_R1.fastq.gz ./${i}_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ./${i}_ \
     --outFilterMultimapNmax 10 \
     --outMultimapperOrder Random \
     --outFilterMismatchNmax 2 \
     --outSAMattributes All \
     --twopassMode Basic --twopass1readsN -1 \
     --quantMode TranscriptomeSAM GeneCounts
