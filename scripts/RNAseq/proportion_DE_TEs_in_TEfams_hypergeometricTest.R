#!/applications/R/R-3.4.0/bin/Rscript

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Perform hypergeometric test to determine whether
#  a significant proportion of elements belonging to each transposon
#  superfamily are differentially expressed (DE) in kyp suvh5 suvh6.
#  P-value is the probability of drawing >= length(DE_TEfamIDs) [x] TEs
#  in a sample size of length(DE_TEIDs) [k] from a total TE set consisting of
#  length(TAIR10_TEfamIDs) [m] + (length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs)) [n]

# Usage 
# /applications/R/R-3.4.0/bin/Rscript ./proportion_DE_TEs_in_TEfams_hypergeometricTest.R ./FDR0.01/res_kssVwt_0.01_lfcShrink_Chr_upRegTEIDs.txt upReg "up-regulated" 100000 100,000

library(methods)

args <- commandArgs(trailingOnly = TRUE)
DE_TEIDsFile <- args[1]
reg <- as.character(args[2])
regulated <- as.character(args[3])
samplesNum <- as.numeric(args[4])
samplesChar <- as.character(args[5])

parDir <- paste0(dirname(DE_TEIDsFile), "/hypergeometricTests/")
outDir <- paste0(parDir, reg, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", parDir, " ] || mkdir ", parDir))
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

TAIR10_TEIDs <- as.character(read.table("/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",
                                        header = T)$Transposon_Name) 

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
TEfamNames <- c(DNAfamNames, RNAfamNames)
TEfamNamesPlot <- c("DNA", "Helitron", "Pogo/Tc1/Mariner", "MuDR", "EnSpm", "hAT", "Harbinger",
                    "RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

# DNA TEs
TAIR10_DNAfamIDs <- lapply(seq_along(DNAfamNames), function(x) {
  as.character(read.table(paste0(DNAdir,
                                 "TAIR10_Buisine_TEs_strand_tab_ann_",
                                 DNAfamNames[x],
                                 ".txt"),
                                 header = T)$name)
})
# RNA TEs
TAIR10_RNAfamIDs <- lapply(seq_along(RNAfamNames), function(x) {
  as.character(read.table(paste0(RNAdir,
                                 "TAIR10_Buisine_TEs_strand_tab_ann_",
                                 RNAfamNames[x],
                                 ".txt"),
                                 header = T)$name)
})

TAIR10_TEfamIDs <- c(TAIR10_DNAfamIDs, TAIR10_RNAfamIDs)

DE_TEIDs <- as.character(read.table(DE_TEIDsFile)$x)

# WARNING: assumes input file has 3-letter extension
baseName <- basename(DE_TEIDsFile)
baseName <- substr(baseName, 1, nchar(baseName)-4)

# P-value is the probability of drawing >= length(DE_TEfamIDs) [x] TEs
# in a sample size of length(DE_TEIDs) [k] from a total TE set consisting of
# length(TAIR10_TEfamIDs) [m] + (length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs)) [n]

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by the sum of the probabilities of drawing
# length(DE_TEfamIDs) to length(DE_TEIDs)

set.seed(18401)
lapply(seq_along(TAIR10_TEfamIDs), function(x)
{
  # Obtain intersection of IDs of differentially expressed TEs and
  # IDs of DNA and RNA family TEs
  DE_TEfamIDs <- intersect(DE_TEIDs, TAIR10_TEfamIDs[[x]])
  print(length(DE_TEfamIDs))
  write.table(DE_TEfamIDs,
              file = paste0(outDir,
                            baseName,
                            "_",
                            TEfamNames[x],
                            ".txt"))
  # Calculate proportion of differentially expressed TEs that are DNA or RNA
  # family elements
  proportion_DE_TEfam <- length(DE_TEfamIDs)/length(DE_TEIDs)

  # Calculate the P-values for over-representation and under-representation
  # of TE family elements among differentially expressed TES
  # Over-representation:
  Pval_overrep <- sum(dhyper(x = length(DE_TEfamIDs):length(DE_TEIDs),
                             m = length(TAIR10_TEfamIDs[[x]]),
                             n = length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs[[x]]),
                             k = length(DE_TEIDs)))
  print(Pval_overrep)
    # Or by 1 minus the sum of the probabilities of drawing 0:(length(DE_TEfamIDs)-1)
  print(1 - sum(dhyper(x = 0:(length(DE_TEfamIDs)-1),
                       m = length(TAIR10_TEfamIDs[[x]]),
                       n = length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs[[x]]),
                       k = length(DE_TEIDs))))

  # Under-representation
  Pval_underrep <- phyper(q = length(DE_TEfamIDs),
                          m = length(TAIR10_TEfamIDs[[x]]),
                          n = length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs[[x]]),
                          k = length(DE_TEIDs))
  print(Pval_underrep)

  # Sample without replacement
  hgDist <- rhyper(nn = samplesNum,
                   m = length(TAIR10_TEfamIDs[[x]]),
                   n = length(TAIR10_TEIDs)-length(TAIR10_TEfamIDs[[x]]),
                   k = length(DE_TEIDs))
  random_proportions_DE_TEfam <- hgDist/length(DE_TEIDs)

  setClass("hypergeomTest",
           representation(Pval_overrep = "numeric",
                          Pval_underrep = "numeric",
                          proportion_DE_TEfam = "numeric",
                          random_proportions_DE_TEfam = "numeric",
                          DE_TEfamIDsLength = "numeric",
                          hypergeometricDistribution = "numeric"))
  hgTestResults <- new("hypergeomTest",
                       Pval_overrep = Pval_overrep,
                       Pval_underrep = Pval_underrep,
                       proportion_DE_TEfam = proportion_DE_TEfam,
                       random_proportions_DE_TEfam = random_proportions_DE_TEfam,
                       DE_TEfamIDsLength = length(DE_TEfamIDs),
                       hypergeometricDistribution = hgDist)
  save(hgTestResults,
       file = paste0(outDir,
                     baseName,
                     "_hypergeometricTestResults_",
                     TEfamNames[x],
                     ".RData"))

  library(plotrix)
  # Generate histogram
  pdf(paste0(plotDir,
             "hist_",
             baseName,
             "_hypergeometricTestResults_",
             TEfamNames[x],
             ".pdf"),
             height = 4, width = 5)
  ## Disable scientific notation (e.g., 0.0001 rather than 1e-04)
  #options(scipen = 100)
  # Calculate max density
  maxDensityPlus <- max(density(hgTestResults@random_proportions_DE_TEfam)$y)*1.2
  if(hgTestResults@proportion_DE_TEfam > mean(hgTestResults@random_proportions_DE_TEfam)) {
    alpha0.05 <- quantile(hgTestResults@random_proportions_DE_TEfam, 0.95)[[1]]
    xlim <- c(pmax(0, min(hgTestResults@random_proportions_DE_TEfam)-.1),
              pmax(hgTestResults@proportion_DE_TEfam+.1, alpha0.05+.1))
    Pval <- Pval_overrep
  } else {
    alpha0.05 <- quantile(hgTestResults@random_proportions_DE_TEfam, 0.05)[[1]]
    xlim <- c(pmax(0, pmin(hgTestResults@proportion_DE_TEfam-.1, alpha0.05-.1)),
              max(hgTestResults@random_proportions_DE_TEfam)+.1)
    Pval <- Pval_underrep
  }
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  hist(hgTestResults@random_proportions_DE_TEfam,
       freq = FALSE,
       col = "grey70",
       border = NA,
       lwd = 2,
       xlim = xlim,
       ylim = c(0,
                maxDensityPlus),
       xlab = paste0("Proportion of ", regulated, " transposable elements \n that are ",
                     TEfamNamesPlot[x], " elements"),
       ylab = "Density",
       main = "",
       cex.lab = 1, cex.axis = 1)
  titleText <- list(bquote(.(baseName)),
                    bquote(italic("P")~" = "~.(Pval)),
                    bquote("Samples (hypergeometric distribution) = "~.(samplesChar)))
  mtext(do.call(expression, titleText), side = 3, line = 2:0, cex = 1)
  lines(density(hgTestResults@random_proportions_DE_TEfam), lwd = 1.5)
  ablineclip(v = mean(hgTestResults@random_proportions_DE_TEfam),
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
  ablineclip(v = hgTestResults@proportion_DE_TEfam,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
  ablineclip(v = alpha0.05,
             y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
  text(x = c(pmax(0.05, min(hgTestResults@random_proportions_DE_TEfam)-.05),
             mean(hgTestResults@random_proportions_DE_TEfam),
             hgTestResults@proportion_DE_TEfam,
             alpha0.05),
       y = c(maxDensityPlus*.95,
             maxDensityPlus,
             maxDensityPlus,
             maxDensityPlus*.95),
       labels = c("Simulated",
                  "Expected",
                  "Observed",
                  expression(alpha~" = 0.05")),
       col = c("grey70",
               "black",
               "forestgreen",
               "red"),
       cex = 0.7)
  box(lwd = 2)
  dev.off()

})
