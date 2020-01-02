#!/applications/R/R-3.3.2/bin/Rscript

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Generate DNA methylation proportion matrices for transposable elements
#  that are upregulated in kyp suvh5 suvh6, for those that are
#  not differentially expressed in kyp suvh5 suvh6, and for random loci.
#  Each row in a matrix represents a feature (a TE or a random locus) and
#  each column represents a window within or flanking the feature

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "/applications/R/R-3.3.2/bin/Rscript ./TE_DNA_methylation_profiles.R 2000 2kb 20 /home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed CGmeth 0.01"

library(EnrichedHeatmap)
library(genomation)

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
covDatPath <- as.character(args[4])
libName <- as.character(args[5])
FDRchar <- as.character(args[6])

# Function to create DNA methylation matrices for feature loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
DNAmethMatrix <- function(signal,
                          feature,
                          ranLoc,
                          featureSize,
                          flankSize,
                          winSize,
                          DNAmethOutDF,
                          DNAmethOutDFcolMeans) {
  # feature loci
  set.seed(2840)
  feature_smoothed <- normalizeToMatrix(signal = signal,
                                        target = feature,
                                        value_column = "coverage",
                                        extend = flankSize,
                                        mean_mode = "absolute",
                                        w = winSize,
                                        empty_value = NA,
                                        smooth = TRUE,
                                        include_target = TRUE,
                                        target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("feature_smoothed")
  print(feature_smoothed)
  print("feature_smoothed rows = ")
  print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  feature_smoothed_failed_rows <- attr(feature_smoothed, "failed_rows")
  print("feature_smoothed failed rows = ")
  print(length(feature_smoothed_failed_rows))
  ## Below code chunk should be commented out to retain complete matrices for
  ## heatmap row ordering
  #if(is.null(feature_smoothed_failed_rows) == FALSE) {
  #  feature_smoothed <- feature_smoothed[-feature_smoothed_failed_rows,]
  #}
  #print("feature_smoothed rows less failed rows = ")
  #print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  ##
  feature_smoothed_DF <- data.frame(feature_smoothed)
  feature_smoothed_DF_colMeans <- as.vector(colMeans(feature_smoothed_DF,
                                                     na.rm = T))
  write.table(feature_smoothed_DF,
              file = DNAmethOutDF[[1]],
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(feature_smoothed_DF_colMeans,
              file = DNAmethOutDFcolMeans[[1]],
              quote = F, sep = "\t", row.names = F, col.names = T)

  # random loci
  set.seed(8472)
  ranLoc_smoothed <- normalizeToMatrix(signal = signal,
                                       target = ranLoc,
                                       value_column = "coverage",
                                       extend = flankSize,
                                       mean_mode = "absolute",
                                       w = winSize,
                                       empty_value = NA,
                                       smooth = TRUE,
                                       include_target = TRUE,
                                       target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("ranLoc_smoothed")
  print(ranLoc_smoothed)
  print("ranLoc_smoothed rows = ")
  print(length(ranLoc_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  ranLoc_smoothed_failed_rows <- attr(ranLoc_smoothed, "failed_rows")
  print("ranLoc_smoothed failed rows = ")
  print(length(ranLoc_smoothed_failed_rows))
  ## Below code chunk should be commented out to retain complete matrices for
  ## heatmap row ordering
  #if(is.null(ranLoc_smoothed_failed_rows) == FALSE) {
  #  ranLoc_smoothed <- ranLoc_smoothed[-ranLoc_smoothed_failed_rows,]
  #}
  #print("ranLoc_smoothed rows less failed rows = ")
  #print(length(ranLoc_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
  ##
  ranLoc_smoothed_DF <- data.frame(ranLoc_smoothed)
  ranLoc_smoothed_DF_colMeans <- as.vector(colMeans(ranLoc_smoothed_DF,
                                                    na.rm = T))
  write.table(ranLoc_smoothed_DF,
              file = DNAmethOutDF[[2]],
              quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(ranLoc_smoothed_DF_colMeans,
              file = DNAmethOutDFcolMeans[[2]],
              quote = F, sep = "\t", row.names = F, col.names = T)
}

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
covGR <- readGeneric(libPath, meta.col = list(coverage = 4))
seqlevels(covGR) <- sub("chr", "Chr", seqlevels(covGR))
assign(paste0(libName), covGR)

inDir <- paste0("/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/DESeq2_TEs_analysis_01/FDR",
                FDRchar, "/")
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

# Convert feature coordinates to GRanges object
TEs <- read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                  header = T)
colnames(TEs) <- c("chr", "start", "end", "strand", "ID", "family", "superfamily")
TEs <- rbind(
  data.frame(TEs[TEs$superfamily == "DNA/En-Spm",],
             superfamName = "EnSpm"),
  data.frame(TEs[TEs$superfamily == "DNA/Harbinger",],
             superfamName = "Harbinger"),
  data.frame(TEs[TEs$superfamily == "DNA/HAT",],
             superfamName = "hAT"),
  data.frame(TEs[TEs$superfamily == "RC/Helitron",],
             superfamName = "Helitron"),
  data.frame(TEs[TEs$superfamily == "DNA/MuDR",],
             superfamName = "MuDR"),
  data.frame(TEs[TEs$superfamily == "DNA/Pogo"
               | TEs$superfamily == "DNA/Tc1"
               | TEs$superfamily == "DNA/Mariner",],
             superfamName = "Pogo_Tc1_Mariner"),
  data.frame(TEs[TEs$superfamily == "DNA",],
             superfamName = "Unclassified_DNA"),
  data.frame(TEs[TEs$superfamily == "LTR/Copia",],
             superfamName = "Copia_LTR"),
  data.frame(TEs[TEs$superfamily == "LTR/Gypsy",],
             superfamName = "Gypsy_LTR"),
  data.frame(TEs[TEs$superfamily == "LINE/L1",],
             superfamName = "LINE1"),
  data.frame(TEs[TEs$superfamily == "LINE?",],
             superfamName = "Putative_LINE"),
  data.frame(TEs[TEs$superfamily == "SINE"
               | TEs$superfamily == "RathE1_cons"
               | TEs$superfamily == "RathE2_cons"
               | TEs$superfamily == "RathE3_cons",],
             superfamName = "SINE"),
  data.frame(TEs[TEs$superfamily == "Unassigned",],
             superfamName = "Unclassified")
)

# Load IDs of differentially expressed TEs
targets <- data.frame(rownames(read.table(paste0(inDir,
                                                "res_kssVwt_",
                                                FDRchar,
                                                "_lfcShrink_Chr_upRegSortedDF_TEs.txt"))))
colnames(targets) <- "ID"
# Extract features whose names match those in "targets"
DEtargets <- TEs[TEs$ID %in% targets$ID,]
superfamNames <- as.character(unique(DEtargets$superfamName))

# Calculate average coverage profile for DE TEs in each superfamily
for(x in seq_along(superfamNames)) {
  print(superfamNames[x])
  superfamDEtargets <- DEtargets[DEtargets$superfamName == superfamNames[x],]
  superfamTEs <- TEs[TEs$superfamName == superfamNames[x],]
  targetsGR <- GRanges(seqnames = superfamDEtargets$chr,
                       ranges = IRanges(start = superfamDEtargets$start,
                                        end = superfamDEtargets$end),
                       strand = superfamDEtargets$strand)
  print(length(targetsGR))

  # Generate GRanges object containing random loci of same number
  # and size distribution as targetsGR
  # Define function to randomly select start coordinates,
  # with the same number per chromosome as targets
  ranLocStartSelect <- function(coordinates, n) {
    sample(x = coordinates,
           size = n,
           replace = FALSE)
  }
  # Define seed so that random selections are reproducible
  set.seed(374592)
  ranLocGR <- GRanges()
  for(i in 1:length(chrs)) {
    targetsGRchr <- targetsGR[seqnames(targetsGR) == chrs[i]]
    if(length(targetsGRchr) > 0) {
      ranLocStartchr <- ranLocStartSelect(coordinates = c((flankSize+max(width(targetsGRchr))+1) :
                                                          (chrLens[i]-max(width(targetsGRchr))-flankSize)),
                                          n = length(targetsGRchr))
      ranLocGRchr <- GRanges(seqnames = chrs[i],
                             ranges = IRanges(start = ranLocStartchr,
                                              width = width(targetsGRchr)),
                             strand = strand(targetsGRchr))
      ranLocGR <- append(ranLocGR, ranLocGRchr)
   }
  }
  stopifnot(identical(sort(width(targetsGR)), sort(width(ranLocGR))))

  # Extract features whose names do not match those of
  # up-regulated features (defined with FDR < 0.1)
  # OR do not match those of down-regulated features (defined with FDR < 0.1)
  # (i.e., features that are not differentially expressed)
  upReg <- data.frame(rownames(read.table(paste0("/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/DESeq2_TEs_analysis_01/FDR0.1/",
                                                 "res_kssVwt_",
                                                 "0.1",
                                                 "_lfcShrink_Chr_upRegSortedDF_TEs.txt"))))
  colnames(upReg) <- "ID"
  downReg <- data.frame(rownames(read.table(paste0("/home/ajt200/analysis/170918_Chris_RNAseq_Col_kss/fastq_pooled/DESeq2_TEs_analysis_01/FDR0.1/",
                                                   "res_kssVwt_",
                                                   "0.1",
                                                   "_lfcShrink_Chr_downRegSortedDF_TEs.txt"))))
  colnames(downReg) <- "ID"
  nonDEtargets <- subset(superfamTEs, !(ID %in% upReg$ID) &
                                      !(ID %in% downReg$ID))
  # Convert nonDEtargets to GRanges object
  nonDEtargetsGR <- GRanges(seqnames = nonDEtargets$chr,
                            ranges = IRanges(start = nonDEtargets$start,
                                             end = nonDEtargets$end),
                            strand = nonDEtargets$strand)
  print(length(nonDEtargetsGR))

  # Define function to randomly select ranges,
  # with the same number per chromosome as targets
  ranNonDEtargetsSelect <- function(nonDEtargetsGR, n) {
    sample(x = nonDEtargetsGR,
           size = n,
           replace = FALSE)
  }

  # Apply ranNonDEtargetsSelect() function on a per-chromosome basis
  # and append the selected ranges to a growing GRanges object
  # Use set.seed() so that random selections can be reproduced
  set.seed(374592)
  ranNonDEtargetsGR <- GRanges()
  for(i in 1:length(chrs)) {
    nonDEtargetsGRchr <- nonDEtargetsGR[seqnames(nonDEtargetsGR) == chrs[i]]
    targetsGRchr <- targetsGR[seqnames(targetsGR) == chrs[i]]
    if(length(targetsGRchr) > 0) {
      ranNonDEtargetsGRchr <- ranNonDEtargetsSelect(nonDEtargetsGR = nonDEtargetsGRchr,
                                                    n = length(targetsGRchr))
      ranNonDEtargetsGR <- append(ranNonDEtargetsGR, ranNonDEtargetsGRchr)
    }
  }
  print(length(ranNonDEtargetsGR))

  # Define matrix and column mean coverage outfile (mean profiles)
  outDF <- list(paste0(matDir, libName,
                       "_norm_cov_", superfamNames[x], "_TEs_smoothed_target_and_", flankName, "_flank_dataframe.txt"),
                paste0(matDir, libName,
                       "_norm_cov_", superfamNames[x], "_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe.txt"))
  outDFcolMeans <- list(paste0(matDir, libName,
                               "_norm_cov_", superfamNames[x], "_TEs_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
                        paste0(matDir, libName,
                               "_norm_cov_", superfamNames[x], "_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))

  # Run DNAmethMatrix() function on each coverage GRanges object to obtain matrices
  ## containing normalised coverage values around target and random loci
  DNAmethMatrix(signal = covGR,
                feature = targetsGR,
                ranLoc = ranNonDEtargetsGR,
                featureSize = mean(width(targetsGR)),
                flankSize = flankSize,
                winSize = winSize,
                DNAmethOutDF = outDF,
                DNAmethOutDFcolMeans = outDFcolMeans)
  print(paste0(libName, " around ", superfamNames[x], " TEs profile calculation complete"))
}
