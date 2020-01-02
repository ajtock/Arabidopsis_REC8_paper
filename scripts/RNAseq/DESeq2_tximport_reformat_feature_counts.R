#!/applications/R/R-3.4.0/bin/Rscript

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Convert transposon-level counts into format for use with DESeq2

library(tximport)
print(packageVersion("tximport"))
#[1] ‘1.4.0’

inDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled"
outDir <- "/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/DESeq2_TEs_analysis_01/"

# Read in table of sample IDs that will be used to specify paths to count files
samples <- read.table(file.path(inDir, "/samples_DESeq2_TEs.txt"), header = T)
print(samples)
#  genotype                     directory                      sample
#1       WT salmon_quants_TEs_strandAware  WT_RNAseq_Chris_Rep1_quant
#2       WT salmon_quants_TEs_strandAware  WT_RNAseq_Chris_Rep2_quant
#3      kss salmon_quants_TEs_strandAware kss_RNAseq_Chris_Rep1_quant
#4      kss salmon_quants_TEs_strandAware kss_RNAseq_Chris_Rep2_quant

# Specify paths to count files
files <- file.path(inDir, samples$genotype, samples$directory, samples$sample, "quant.sf")
# Set 1:#_samples
names(files) <- paste0("sample", 1:4)
all(file.exists(files))

# Create a dataframe of transcript IDs and corresponding TE IDs
transID <- read.table(files[1], colClasses = c(NA, rep("NULL", 4)), header = T)
tx2gene <- data.frame(cbind(as.vector(transID[,1]), as.vector(transID[,1])))
colnames(tx2gene) <- c("TXNAME", "TXNAME2")

# Import transcript-level counts, summarised at transcript level with "txOut = TRUE"
library(readr)
txi <- tximport(files, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
print(names(txi))
#[1] "abundance"           "counts"              "length"             
#[4] "countsFromAbundance"

print(head(txi$counts))
#             sample1   sample2    sample3    sample4
#AT1TE52125 64.816153 84.826588 113.638349 110.731785
#AT1TE42735  2.822714  5.088586  16.037173   6.081396
#AT1TE36140 53.651001 70.493246  43.819152  51.017245
#AT1TE21850 18.000000  2.173260   4.688952   5.000000
#AT1TE95105  2.000000  1.000000   1.083565   2.000000
#AT1TE18375 42.042256 38.338420  35.115370  36.543303

save(txi,
     file = paste0(outDir, "tximport.RData"))

