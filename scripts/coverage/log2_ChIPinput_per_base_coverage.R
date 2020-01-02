#!/applications/R/R-3.5.0/bin/Rscript

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Calculate log2 ratio of ChIP to input coverage
#  for each based and write to 1-based BED-like file suitable
#  for use with  the Bioconductor package EnrichedHeatmap

# Usage via Condor submission system on hydrogen node7:
# csmit -m 20G -c 1 "/applications/R/R-3.5.0/bin/Rscript ./log2_ChIPinput_per_base_coverage.R wt_REC8_HA_Rep2_ChIP wt_REC8_Myc_Rep1_input"

# Note: wt_REC8_Myc_Rep1_input is used as a common input for
# log2 transformation of all wild-type REC8 libraries, as this
# input library was sequenced at the greatest depth and therefore
# had the most uniform coverage distribution across the genome

args <- commandArgs(trailingOnly = TRUE)
ChIPname <- args[1]
inputname <- args[2]
outname <- paste0("log2_", ChIPname, "_", inputname)

ChIP <- read.table(paste0(ChIPname,
                          "_norm_allchrs_coverage_coord_tab.bed"))
print(head(ChIP))
input <- read.table(paste0(inputname,
                           "_norm_allchrs_coverage_coord_tab.bed"))
print(head(input))
ChIP_input_offset <- cbind(ChIP, ChIP$V4+1, input$V4+1)  
ChIP_input_offset <- ChIP_input_offset[,-4]
colnames(ChIP_input_offset) <- c("V1", "V2", "V3", "V4", "V5")
print(head(ChIP_input_offset))
norm <- log2(ChIP_input_offset$V4/ChIP_input_offset$V5)
ChIPnormInput <- cbind(ChIP_input_offset, norm)
log2ChIPinput <- ChIPnormInput[,-4:-5] 
write.table(log2ChIPinput,
            file = paste0("./",
                          outname,
                          "_norm_allchrs_coverage_coord_tab.bed"),
            sep = "\t", quote = F,
            row.names = F, col.names = F)
rm(ChIPnormInput)
rm(log2ChIPinput)
gc()

