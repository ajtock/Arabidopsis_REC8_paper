#!/bin/bash

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Concatenate per-chromosome bedgraphs generated by
#  library_size_normalised_per_base_coverage.R and convert
#  into BED-like format (but with 1-based coordinates) suitable for
#  use with the Bioconductor package EnrichedHeatmap

# Example usage via condor submission system on hydrogen node7:
# csmit -m 10G -c 1 "bash ./concatenate_chromosome_bedgraphs.sh wt_REC8_HA_Rep2_ChIP"

prefix=$1

# Concatenate chromosome bedgraph files to create genome bedgraph file (0-based)
cat $(find ./ -name ${prefix}"_norm_Chr*_coverage.bedgraph" | sort -V) > ${prefix}_libsize_norm_cov.bedGraph
sed -i '1i track type=bedGraph' ${prefix}_libsize_norm_cov.bedGraph
# Convert 0-based genome bedgraph file to 1-based bed-like file containing coordinates and library-normalised coverage
awk 'BEGIN {OFS="\t"}; {print $1, $3, $3, $4}' ${prefix}_libsize_norm_cov.bedGraph > ${prefix}_norm_allchrs_coverage_coord_tab_tmp.bed
tail -n +2 ${prefix}_norm_allchrs_coverage_coord_tab_tmp.bed > ${prefix}_norm_allchrs_coverage_coord_tab.bed
rm ${prefix}_norm_allchrs_coverage_coord_tab_tmp.bed
