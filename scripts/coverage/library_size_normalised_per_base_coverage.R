#!/applications/R/R-3.3.2/bin/Rscript

# Author: Xiaohui Zhao (xz289@cam.ac.uk)
# Description: Calculate per-base coverage normalised by
#  total genome-wide coverage

# Usage via Condor submission system on hydrogen node7:
# csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./library_size_normalised_per_base_coverage.R wt_REC8_HA_Rep2_ChIP"

args <- commandArgs(trailingOnly = TRUE)
libName1 <- args[1]
libFile <- paste0(libName1, "_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam")

library(Rsamtools)
library(GenomicAlignments)

# Define function to extract coverage from bowtie2-generated BAM file
# derived from paired-end data
chrdat_pair_function <- function(index, libfiles, covfiles){
  un1 <- libfiles[index]
  reads1 <- readGAlignmentPairs(un1)
  cov1 <- coverage(reads1)
  chrs <- list()
  chrs[[1]] <- as.integer(cov1$'1'); chrs[[2]] <- as.integer(cov1$'2')
  chrs[[3]] <- as.integer(cov1$'3'); chrs[[4]] <- as.integer(cov1$'4')
  chrs[[5]] <- as.integer(cov1$'5')
  for(a in 1:5){
    out.file <- covfiles[[index]][a]
    write.table(chrs[[a]], file = out.file, row.names = F, col.names = F)
    gc()
  }
}

chrIndex <- c(1:5)
covFiles <- list(paste0("./", libName1, "_chr", chrIndex, "_unique_both_coverage.txt"))

for(i in 1:length(libName1)) {
  print(i)
  index <- i
  test <- chrdat_pair_function(index = i,
                               libfiles = libFile,
                               covfiles = covFiles)
  gc()
}

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrLens <- c(30427671,19698289,23459830,18585056,26975502)

j <- 1
print(libName1)
outFiles <- paste0("./", libName1, "_norm_Chr", c(1:5), "_coverage.bedgraph")
covFilesVec <- covFiles[[j]]
myFiles1 <- read.table(covFilesVec[1], header=F)
myFiles2 <- read.table(covFilesVec[2], header=F)
myFiles3 <- read.table(covFilesVec[3], header=F)
myFiles4 <- read.table(covFilesVec[4], header=F)
myFiles5 <- read.table(covFilesVec[5], header=F)
totalSum <- colSums(myFiles1)+colSums(myFiles2)+colSums(myFiles3)+colSums(myFiles4)+colSums(myFiles5)
##
myFiles1.new1 <- myFiles1*1e9/totalSum
myFiles1.new <- cbind(chrs[1], c(0:(chrLens[1]-1)), c(1:chrLens[1]), myFiles1.new1)
write.table(myFiles1.new, file=outFiles[1], row.names=F, col.names=F, quote=F)
##
myFiles2.new1 <- myFiles2*1e9/totalSum
myFiles2.new <- cbind(chrs[2], c(0:(chrLens[2]-1)), c(1:chrLens[2]), myFiles2.new1)
write.table(myFiles2.new, file=outFiles[2], row.names=F, col.names=F, quote=F)
##
myFiles3.new1 <- myFiles3*1e9/totalSum
myFiles3.new <- cbind(chrs[3], c(0:(chrLens[3]-1)), c(1:chrLens[3]), myFiles3.new1)
write.table(myFiles3.new, file=outFiles[3], row.names=F, col.names=F, quote=F)
##
myFiles4.new1 <- myFiles4*1e9/totalSum
myFiles4.new <- cbind(chrs[4], c(0:(chrLens[4]-1)), c(1:chrLens[4]), myFiles4.new1)
write.table(myFiles4.new, file=outFiles[4], row.names=F, col.names=F, quote=F)
##
myFiles5.new1 <- myFiles5*1e9/totalSum
myFiles5.new <- cbind(chrs[5], c(0:(chrLens[5]-1)), c(1:chrLens[5]), myFiles5.new1)
write.table(myFiles5.new, file=outFiles[5], row.names=F, col.names=F, quote=F)

