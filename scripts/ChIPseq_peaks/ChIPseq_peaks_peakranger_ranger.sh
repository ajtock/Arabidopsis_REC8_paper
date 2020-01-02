#!/bin/bash

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Use ranger tool within PeakRanger (version 1.18)
#  to call peaks in a ChIP-seq library
# PeakRanger manual: http://ranger.sourceforge.net/manual1.18.html

# Usage via Condor submission system on hydrogen node7:
# csmit -m 10G -c 8 "bash ./ChIPseq_peaks_peakranger_ranger.sh '/home/ajt200/analysis/REC8_pooled' '/home/ajt200/analysis/REC8_pooled' wt_REC8_HA_Rep2_ChIP wt_REC8_Myc_Rep1_input 0.001 0.01 200 8"

# Note: wt_REC8_Myc_Rep1_input is used as a common input for
# log2 transformation of all wild-type REC8 libraries, as this
# input library was sequenced at the greatest depth and therefore
# had the most uniform coverage distribution across the genome

ChIP_bamDir=$1
input_bamDir=$2
ChIP_prefix=$3
input_prefix=$4
pval=$5
qval=$6
ext_length=$7
threads=$8

/home/ajt200/tools/PeakRanger-1.18/bin/peakranger ranger \
  --data ${ChIP_bamDir}/${ChIP_prefix}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
  --control ${input_bamDir}/${input_prefix}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
  --format bam \
  --output ${ChIP_prefix}"_rangerPeaks" \
  --pval ${pval} \
  --FDR ${qval} \
  --ext_length ${ext_length} \
  --pad \
  --thread ${threads} \
  --verbose
