#!/applications/R/R-3.4.0/bin/Rscript

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Generate coverage matrices for peaks and random loci.
#  Each row in a matrix represents a feature (a peak or a random locus) and
#  each column represents a window within or flanking the feature

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./peak_coverage_profiles.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab.bed REC8_HA_Rep2"

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
covDatPath <- as.character(args[4])
libName <- as.character(args[5])

# Function to create coverage matrices for feature loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
covMatrix <- function(signal,
                      feature,
                      ranLoc,
                      featureSize,
                      flankSize,
                      winSize,
                      outDF,
                      outDFcolMeans) {
  # feature loci
  set.seed(2840)
  feature_smoothed <- normalizeToMatrix(signal = signal,
                                        target = feature,
                                        value_column = "coverage",
                                        extend = flankSize,
                                        mean_mode = "w0",
                                        w = winSize,
                                        background = 0,
                                        smooth = TRUE,
                                        include_target = TRUE,
                                        target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("feature_smoothed")
  print(feature_smoothed)
  print("feature_smoothed rows = ")
  print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  feature_smoothed_DF <- data.frame(feature_smoothed)
  feature_smoothed_DF_colMeans <- as.vector(colMeans(feature_smoothed_DF,
                                                     na.rm = T))
  write.table(feature_smoothed_DF,
              file = outDF[[1]],
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(feature_smoothed_DF_colMeans,
              file = outDFcolMeans[[1]],
              quote = F, sep = "\t", row.names = F, col.names = T)

  # random loci
  set.seed(8472)
  ranLoc_smoothed <- normalizeToMatrix(signal = signal,
                                       target = ranLoc,
                                       value_column = "coverage",
                                       extend = flankSize,
                                       mean_mode = "w0",
                                       w = winSize,
                                       background = 0,
                                       smooth = TRUE,
                                       include_target = TRUE,
                                       target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("ranLoc_smoothed")
  print(ranLoc_smoothed)
  print("ranLoc_smoothed rows = ")
  print(length(ranLoc_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  ranLoc_smoothed_DF <- data.frame(ranLoc_smoothed)
  ranLoc_smoothed_DF_colMeans <- as.vector(colMeans(ranLoc_smoothed_DF,
                                                    na.rm = T))
  write.table(ranLoc_smoothed_DF,
              file = outDF[[2]],
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(ranLoc_smoothed_DF_colMeans,
              file = outDFcolMeans[[2]],
              quote = F, sep = "\t", row.names = F, col.names = T)
}

matDir <- "./matrices/"
plotDir <- "./plots/"
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Chromosome definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(rep(1, 5))
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
mask <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))

# Import peaks as GRanges object
load(paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/",
            "REC8_HA_Rep2_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"))
peaksGR <- rangerPeaksGR_arm_mergedOverlaps
strand(peaksGR) <- "*"
peaksGR <- sortSeqlevels(peaksGR)
peaksGR <- sort(peaksGR)
print("***********peaks***********")
print(peaksGR)
# Generate GRanges object containing random loci of same number
# and size distribution as peaksGR
set.seed(374592)
ranLocGR <- randomizeRegions(peaksGR,
                             genome = genome,
                             mask = mask,
                             per.chromosome = TRUE,
                             allow.overlaps = TRUE)
# Confirm ranLocGR does not contain loci in masked regions
stopifnot(sum(countOverlaps(mask, ranLocGR)) == 0)

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
covGR <- readGeneric(libPath, meta.col = list(coverage = 4))
assign(paste0(libName), covGR)

# Define matrix and column mean coverage outfile (mean profiles)
outDF <- list(paste0(matDir, libName,
                     "_norm_cov_feature_smoothed_target_and_", flankName, "_flank_dataframe.txt"),
              paste0(matDir, libName,
                     "_norm_cov_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, libName,
                             "_norm_cov_feature_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, libName,
                             "_norm_cov_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
covMatrix(signal = covGR,
          feature = peaksGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(peaksGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(libName, " profile calculation complete"))


