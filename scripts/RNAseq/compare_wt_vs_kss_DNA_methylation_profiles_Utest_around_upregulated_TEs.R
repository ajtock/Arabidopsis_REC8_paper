#!/applications/R/R-3.5.0/bin/Rscript

# Author: Andrew Tock (ajt200@cam.ac.uk)
# Description: Use Mann-Whitney-Wilcoxon (U) tests to determine whether significant
#  differences exist between genotypes with regard to DNA methylation proportion values
#  within windows along upregulated transposable elements,
#  non-differentially expressed transposable elements, and random loci.
#  Plot genotype mean DNA methylation proportion profiles with 95% confidence intervals as ribbons,
#  and with vertical ticks along the x-axis indicating windows in which
#  significant differences exist between genotypes

# Usage:
# /applications/R/R-3.5.0/bin/Rscript compare_wt_vs_kss_DNA_methylation_profiles_Utest_around_upregulated_TEs.R Upregulated_TEs_kss 'Upregulated' 2000 2kb '2 kb' 20 20bp 'mCHG,kss_mCHG' 'wt mCHG,kss mCHG' 'orange,orange4' '0.02,0.50' 0.1

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
featureNamePlot <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
libNames <- unlist(strsplit(args[8],
                            split = ","))
libNamesPlot <- unlist(strsplit(args[9],
                                split = ","))
colours <- unlist(strsplit(args[10],
                           split = ","))
legendPos <- as.numeric(unlist(strsplit(args[11],
                                        split = ",")))
FDR <- as.numeric(args[12])
FDRname <- paste0("FDR", as.character(FDR))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

# Define and create directories and subdirectories
matDir <- "./matrices/"
ranMatDir <- "./ranLoc/matrices/"
outDir <- paste0("./ranLoc/Utests_", FDRname, "/")
plotDir <- paste0(outDir, "/plots_sig_wins/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

superfamNames <- system(paste0("ls ", matDir, "REC8_HA_Rep2*TEs*dataframe.txt"),
                        intern = T)
superfamNames <- sub(pattern = "./matrices/REC8_HA_Rep2_norm_cov_(\\w+)_TEs\\w+.txt",
                     replacement = "\\1",
                     x = superfamNames)

for(y in seq_along(superfamNames)) {
  # Load feature coverage matrix for each ChIP-seq dataset
  featureMats <- mclapply(seq_along(libNames), function(x) {
    as.matrix(read.table(paste0(matDir,
                                libNames[x],
                                "_norm_cov_",
                                superfamNames[y],
                                "_TEs_smoothed_target_and_",
                                flankName, "_flank_dataframe.txt"),
                         header = T))
  }, mc.cores = length(libNames))
  # Convert into data.frame and run U-tests
  featureMatsDF <- lapply(seq_along(featureMats), function(x) {
    as.data.frame(featureMats[[x]])
  })
  featureMatsDF_Utests <- lapply(seq_along(featureMatsDF[[1]]), function(x) {
    wilcox.test(x = featureMatsDF[[1]][[x]],
                y = featureMatsDF[[2]][[x]],
                alternative = "two.sided")
  })
  # "cannot compute exact p-value with ties"
  featureMatsDF_UtestsPvals <- sapply(seq_along(featureMatsDF_Utests), function(x) {
    featureMatsDF_Utests[[x]]$p.val
  })
  # Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
  featureMatsDF_UtestsAdjPvals <- p.adjust(p = featureMatsDF_UtestsPvals,
                                           method = "BH")
  featureMatsDF_UtestsSigWins <- which(-log10(featureMatsDF_UtestsAdjPvals) > -log10(FDR))

  # Load random non-DE TEs coverage matrix for each ChIP-seq dataset
  ranLocMats <- mclapply(seq_along(libNames), function(x) {
    as.matrix(read.table(paste0(matDir,
                                libNames[x],
                                "_norm_cov_",
                                superfamNames[y],
                                "_ranLoc_smoothed_target_and_",
                                flankName, "_flank_dataframe.txt"),
                         header = T))
  }, mc.cores = length(libNames))
  # Convert into data.frame and run U-tests
  ranLocMatsDF <- lapply(seq_along(ranLocMats), function(x) {
    as.data.frame(ranLocMats[[x]])
  })
  ranLocMatsDF_Utests <- lapply(seq_along(ranLocMatsDF[[1]]), function(x) {
    wilcox.test(x = ranLocMatsDF[[1]][[x]],
                y = ranLocMatsDF[[2]][[x]],
                alternative = "two.sided")
  })
  # "cannot compute exact p-value with ties"
  ranLocMatsDF_UtestsPvals <- sapply(seq_along(ranLocMatsDF_Utests), function(x) {
    ranLocMatsDF_Utests[[x]]$p.val
  })
  # Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
  ranLocMatsDF_UtestsAdjPvals <- p.adjust(p = ranLocMatsDF_UtestsPvals,
                                          method = "BH")
  ranLocMatsDF_UtestsSigWins <- which(-log10(ranLocMatsDF_UtestsAdjPvals) > -log10(FDR))

  # Load ranLoc coverage matrix for each ChIP-seq dataset
  ranLoc2Mats <- mclapply(seq_along(libNames), function(x) {
    as.matrix(read.table(paste0(ranMatDir,
                                libNames[x],
                                "_norm_cov_",
                                superfamNames[y],
                                "_ranLoc_smoothed_target_and_",
                                flankName, "_flank_dataframe.txt"),
                         header = T))
  }, mc.cores = length(libNames))
  # Convert into data.frame and run U-tests
  ranLoc2MatsDF <- lapply(seq_along(ranLoc2Mats), function(x) {
    as.data.frame(ranLoc2Mats[[x]])
  })
  ranLoc2MatsDF_Utests <- lapply(seq_along(ranLoc2MatsDF[[1]]), function(x) {
    wilcox.test(x = ranLoc2MatsDF[[1]][[x]],
                y = ranLoc2MatsDF[[2]][[x]],
                alternative = "two.sided")
  })
  # "cannot compute exact p-value with ties"
  ranLoc2MatsDF_UtestsPvals <- sapply(seq_along(ranLoc2MatsDF_Utests), function(x) {
    ranLoc2MatsDF_Utests[[x]]$p.val
  })
  # Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
  ranLoc2MatsDF_UtestsAdjPvals <- p.adjust(p = ranLoc2MatsDF_UtestsPvals,
                                           method = "BH")
  ranLoc2MatsDF_UtestsSigWins <- which(-log10(ranLoc2MatsDF_UtestsAdjPvals) > -log10(FDR))

  ## feature
  # Transpose matrix and convert to dataframe
  # in which first column is window name
  wideDFfeature_list <- mclapply(seq_along(featureMats), function(x) {
    data.frame(window = colnames(featureMats[[x]]),
               t(featureMats[[x]]))
  }, mc.cores = length(featureMats))
  
  # Convert into tidy data.frame (long format)
  tidyDFfeature_list  <- mclapply(seq_along(wideDFfeature_list), function(x) {
    gather(data  = wideDFfeature_list[[x]],
           key   = feature,
           value = coverage,
           -window)
  }, mc.cores = length(wideDFfeature_list))
  
  # Order levels of factor "window" so that sequential levels
  # correspond to sequential windows
  for(x in seq_along(tidyDFfeature_list)) {
    tidyDFfeature_list[[x]]$window <- factor(tidyDFfeature_list[[x]]$window,
                                             levels = as.character(wideDFfeature_list[[x]]$window))
  }
  
  # Create summary data.frame in which each row corresponds to a window (Column 1),
  # Column2 is the number of coverage values (features) per window,
  # Column3 is the mean of coverage values per window,
  # Column4 is the standard deviation of coverage values per window,
  # Column5 is the standard error of the mean of coverage values per window,
  # Column6 is the lower bound of the 95% confidence interval, and
  # Column7 is the upper bound of the 95% confidence interval
  summaryDFfeature_list  <- mclapply(seq_along(tidyDFfeature_list), function(x) {
    data.frame(window = as.character(wideDFfeature_list[[x]]$window),
               n      = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                               INDEX = tidyDFfeature_list[[x]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                               INDEX = tidyDFfeature_list[[x]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                               INDEX = tidyDFfeature_list[[x]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  }, mc.cores = length(tidyDFfeature_list))
  
  for(x in seq_along(summaryDFfeature_list)) {
    summaryDFfeature_list[[x]]$window <- factor(summaryDFfeature_list[[x]]$window,
                                                levels = as.character(wideDFfeature_list[[x]]$window))
    summaryDFfeature_list[[x]]$winNo <- factor(1:dim(summaryDFfeature_list[[x]])[1])
    summaryDFfeature_list[[x]]$sem <- summaryDFfeature_list[[x]]$sd/sqrt(summaryDFfeature_list[[x]]$n-1)
    summaryDFfeature_list[[x]]$CI_lower <- summaryDFfeature_list[[x]]$mean -
      qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
    summaryDFfeature_list[[x]]$CI_upper <- summaryDFfeature_list[[x]]$mean +
      qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
  }
  
  names(summaryDFfeature_list) <- libNamesPlot
  
  # Convert list summaryDFfeature_list into a single data.frame for plotting
  summaryDFfeature <- bind_rows(summaryDFfeature_list, .id = "libName")
  summaryDFfeature$libName <- factor(summaryDFfeature$libName,
                                     levels = names(summaryDFfeature_list))
  
  
  ## ranLoc
  # Transpose matrix and convert to dataframe
  # in which first column is window name
  wideDFranLoc_list <- mclapply(seq_along(ranLocMats), function(x) {
    data.frame(window = colnames(ranLocMats[[x]]),
               t(ranLocMats[[x]]))
  }, mc.cores = length(ranLocMats))
  
  # Convert into tidy data.frame (long format)
  tidyDFranLoc_list  <- mclapply(seq_along(wideDFranLoc_list), function(x) {
    gather(data  = wideDFranLoc_list[[x]],
           key   = ranLoc,
           value = coverage,
           -window)
  }, mc.cores = length(wideDFranLoc_list))
  
  # Order levels of factor "window" so that sequential levels
  # correspond to sequential windows
  for(x in seq_along(tidyDFranLoc_list)) {
    tidyDFranLoc_list[[x]]$window <- factor(tidyDFranLoc_list[[x]]$window,
                                            levels = as.character(wideDFranLoc_list[[x]]$window))
  }
  
  # Create summary data.frame in which each row corresponds to a window (Column 1),
  # Column2 is the number of coverage values (ranLocs) per window,
  # Column3 is the mean of coverage values per window,
  # Column4 is the standard deviation of coverage values per window,
  # Column5 is the standard error of the mean of coverage values per window,
  # Column6 is the lower bound of the 95% confidence interval, and
  # Column7 is the upper bound of the 95% confidence interval
  summaryDFranLoc_list  <- mclapply(seq_along(tidyDFranLoc_list), function(x) {
    data.frame(window = as.character(wideDFranLoc_list[[x]]$window),
               n      = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                               INDEX = tidyDFranLoc_list[[x]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                               INDEX = tidyDFranLoc_list[[x]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                               INDEX = tidyDFranLoc_list[[x]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  }, mc.cores = length(tidyDFranLoc_list))
  
  for(x in seq_along(summaryDFranLoc_list)) {
    summaryDFranLoc_list[[x]]$window <- factor(summaryDFranLoc_list[[x]]$window,
                                               levels = as.character(wideDFranLoc_list[[x]]$window))
    summaryDFranLoc_list[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list[[x]])[1])
    summaryDFranLoc_list[[x]]$sem <- summaryDFranLoc_list[[x]]$sd/sqrt(summaryDFranLoc_list[[x]]$n-1)
    summaryDFranLoc_list[[x]]$CI_lower <- summaryDFranLoc_list[[x]]$mean -
      qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
    summaryDFranLoc_list[[x]]$CI_upper <- summaryDFranLoc_list[[x]]$mean +
      qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
  }
  
  names(summaryDFranLoc_list) <- libNamesPlot
  
  # Convert list summaryDFranLoc_list into a single data.frame for plotting
  summaryDFranLoc <- bind_rows(summaryDFranLoc_list, .id = "libName")
  summaryDFranLoc$libName <- factor(summaryDFranLoc$libName,
                                    levels = names(summaryDFranLoc_list))
  

  ## ranLoc2
  # Transpose matrix and convert to dataframe
  # in which first column is window name
  wideDFranLoc2_list <- mclapply(seq_along(ranLoc2Mats), function(x) {
    data.frame(window = colnames(ranLoc2Mats[[x]]),
               t(ranLoc2Mats[[x]]))
  }, mc.cores = length(ranLoc2Mats))
  
  # Convert into tidy data.frame (long format)
  tidyDFranLoc2_list  <- mclapply(seq_along(wideDFranLoc2_list), function(x) {
    gather(data  = wideDFranLoc2_list[[x]],
           key   = ranLoc2,
           value = coverage,
           -window)
  }, mc.cores = length(wideDFranLoc2_list))
  
  # Order levels of factor "window" so that sequential levels
  # correspond to sequential windows
  for(x in seq_along(tidyDFranLoc2_list)) {
    tidyDFranLoc2_list[[x]]$window <- factor(tidyDFranLoc2_list[[x]]$window,
                                             levels = as.character(wideDFranLoc2_list[[x]]$window))
  }
  
  # Create summary data.frame in which each row corresponds to a window (Column 1),
  # Column2 is the number of coverage values (ranLoc2s) per window,
  # Column3 is the mean of coverage values per window,
  # Column4 is the standard deviation of coverage values per window,
  # Column5 is the standard error of the mean of coverage values per window,
  # Column6 is the lower bound of the 95% confidence interval, and
  # Column7 is the upper bound of the 95% confidence interval
  summaryDFranLoc2_list  <- mclapply(seq_along(tidyDFranLoc2_list), function(x) {
    data.frame(window = as.character(wideDFranLoc2_list[[x]]$window),
               n      = tapply(X     = tidyDFranLoc2_list[[x]]$coverage,
                               INDEX = tidyDFranLoc2_list[[x]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFranLoc2_list[[x]]$coverage,
                               INDEX = tidyDFranLoc2_list[[x]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFranLoc2_list[[x]]$coverage,
                               INDEX = tidyDFranLoc2_list[[x]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  }, mc.cores = length(tidyDFranLoc2_list))
  
  for(x in seq_along(summaryDFranLoc2_list)) {
    summaryDFranLoc2_list[[x]]$window <- factor(summaryDFranLoc2_list[[x]]$window,
                                               levels = as.character(wideDFranLoc2_list[[x]]$window))
    summaryDFranLoc2_list[[x]]$winNo <- factor(1:dim(summaryDFranLoc2_list[[x]])[1])
    summaryDFranLoc2_list[[x]]$sem <- summaryDFranLoc2_list[[x]]$sd/sqrt(summaryDFranLoc2_list[[x]]$n-1)
    summaryDFranLoc2_list[[x]]$CI_lower <- summaryDFranLoc2_list[[x]]$mean -
      qt(0.975, df = summaryDFranLoc2_list[[x]]$n-1)*summaryDFranLoc2_list[[x]]$sem
    summaryDFranLoc2_list[[x]]$CI_upper <- summaryDFranLoc2_list[[x]]$mean +
      qt(0.975, df = summaryDFranLoc2_list[[x]]$n-1)*summaryDFranLoc2_list[[x]]$sem
  }
  
  names(summaryDFranLoc2_list) <- libNamesPlot
  
  # Convert list summaryDFranLoc2_list into a single data.frame for plotting
  summaryDFranLoc2 <- bind_rows(summaryDFranLoc2_list, .id = "libName")
  summaryDFranLoc2$libName <- factor(summaryDFranLoc2$libName,
                                     levels = names(summaryDFranLoc2_list))
  
  featureStartLab <- "Start"
  featureEndLab <- "End"
  
  # Define y-axis limits
  #ymin <- min(c(summaryDFfeature$mean-summaryDFfeature$sem,
  #              summaryDFranLoc$mean-summaryDFranLoc$sem))
  #ymax <- max(c(summaryDFfeature$mean+summaryDFfeature$sem,
  #              summaryDFranLoc$mean+summaryDFranLoc$sem))
  ymin <- min(c(summaryDFfeature$CI_lower,
                summaryDFranLoc$CI_lower,
                summaryDFranLoc2$CI_lower))
  ymax <- max(c(summaryDFfeature$CI_upper,
                summaryDFranLoc$CI_upper,
                summaryDFranLoc2$CI_upper))
  
  # Function for formatting y-axis labels
  # with a given number of decimals
  fmt_decimals <- function(decimals) {
    function(x) format(x, nsmall = decimals, scientific = FALSE)
  }
  
  # Define legend labels
  legendLab1 <- grobTree(textGrob(bquote(.(libNamesPlot[1])),
                                  x = legendPos[1], y = legendPos[2], just = "left",
                                  gp = gpar(col = colours[1], fontsize = 18))) 
  legendLab2 <- grobTree(textGrob(bquote(italic(.(unlist(strsplit(libNamesPlot[2],
                                                                  split = " "))[1])) ~
                                         .(unlist(strsplit(libNamesPlot[2],
                                                           split = " "))[2])),
                                  x = legendPos[1], y = legendPos[2]-0.1, just = "left",
                                  gp = gpar(col = colours[2], fontsize = 18)))
  
  # Plot average coverage profiles with 95% CI ribbon
  ## feature
  ggObjGA <- NULL
  ggObj1 <- NULL
  ggObj1 <- ggplot(data = summaryDFfeature,
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = libName),
                  ) +
    geom_rug(data = as.data.frame(featureMatsDF_UtestsSigWins),
             mapping = aes(x = featureMatsDF_UtestsSigWins),
             sides = "b",
             colour = colours[2],
             inherit.aes = FALSE) +
    geom_line(data = summaryDFfeature,
              mapping = aes(colour = libName),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFfeature,
                #mapping = aes(ymin = mean-sem,
                #              ymax = mean+sem,
                mapping = aes(ymin = CI_lower,
                              ymax = CI_upper,
                              fill = libName),
                alpha = 0.4) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(ymin, ymax),
                       labels = function(x) sprintf("%5.2f", x)) +
    scale_x_discrete(breaks = c(1,
                                (upstream/binSize)+1,
                                (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize),
                                dim(summaryDFfeature_list[[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                featureStartLab,
                                featureEndLab,
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize)),
               linetype = "dashed",
               size = 1) +
    labs(x = "",
         y = "DNA methylation") +
    annotation_custom(legendLab1) +
    annotation_custom(legendLab2) +
    theme_bw() +
    theme(
          axis.ticks = element_line(size = 1.0, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size = 22, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 30, colour = "black"),
          legend.position = "none",
          #legend.text = element_text(size = 10),
          #legend.background = element_rect(fill = "transparent"),
          #legend.key = element_rect(colour = "transparent",
          #                          fill = "transparent"),
          #legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 3.5, colour = "black"),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote(.(featureNamePlot) ~ .(superfamNames[y]) ~
                   "(" * italic("n") ~ "=" ~
                   .(prettyNum(summaryDFfeature$n[1],
                               big.mark = ",", trim = T)) *
                   ")"))

  ## ranLoc
  ggObj2 <- NULL
  ggObj2 <- ggplot(data = summaryDFranLoc,
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = libName),
                  ) +
    geom_rug(data = as.data.frame(ranLocMatsDF_UtestsSigWins),
             mapping = aes(x = ranLocMatsDF_UtestsSigWins),
             sides = "b",
             colour = colours[2],
             inherit.aes = FALSE) +
    geom_line(data = summaryDFranLoc,
              mapping = aes(colour = libName),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFranLoc,
                #mapping = aes(ymin = mean-sem,
                #              ymax = mean+sem,
                mapping = aes(ymin = CI_lower,
                              ymax = CI_upper,
                              fill = libName),
                alpha = 0.4) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(ymin, ymax),
                       labels = function(x) sprintf("%5.2f", x)) +
    scale_x_discrete(breaks = c(1,
                                (upstream/binSize)+1,
                                (dim(summaryDFranLoc_list[[1]])[1])-(downstream/binSize),
                                dim(summaryDFranLoc_list[[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                "Start",
                                "End",
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFranLoc_list[[1]])[1])-(downstream/binSize)),
               linetype = "dashed",
               size = 1) +
    labs(x = "",
         y = "") +
    annotation_custom(legendLab1) +
    annotation_custom(legendLab2) +
    theme_bw() +
    theme(
          axis.ticks = element_line(size = 1.0, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size = 22, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 30, colour = "black"),
          legend.position = "none",
          #legend.text = element_text(size = 10),
          #legend.background = element_rect(fill = "transparent"),
          #legend.key = element_rect(colour = "transparent",
          #                          fill = "transparent"),
          #legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 3.5, colour = "black"),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote("Non-DE" ~ .(superfamNames[y]) ~
                   "(" * italic("n") ~ "=" ~
                   .(prettyNum(summaryDFranLoc$n[1],
                               big.mark = ",", trim = T)) *
                   ")"))
  ## ranLoc2
  ggObj3 <- NULL
  ggObj3 <- ggplot(data = summaryDFranLoc2,
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = libName),
                  ) +
    geom_rug(data = as.data.frame(ranLoc2MatsDF_UtestsSigWins),
             mapping = aes(x = ranLoc2MatsDF_UtestsSigWins),
             sides = "b",
             colour = colours[2],
             inherit.aes = FALSE) +
    geom_line(data = summaryDFranLoc2,
              mapping = aes(colour = libName),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFranLoc2,
                #mapping = aes(ymin = mean-sem,
                #              ymax = mean+sem,
                mapping = aes(ymin = CI_lower,
                              ymax = CI_upper,
                              fill = libName),
                alpha = 0.4) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(ymin, ymax),
                       labels = function(x) sprintf("%5.2f", x)) +
    scale_x_discrete(breaks = c(1,
                                (upstream/binSize)+1,
                                (dim(summaryDFranLoc2_list[[1]])[1])-(downstream/binSize),
                                dim(summaryDFranLoc2_list[[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                "Start",
                                "End",
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFranLoc2_list[[1]])[1])-(downstream/binSize)),
               linetype = "dashed",
               size = 1) +
    labs(x = "",
         y = "") +
    annotation_custom(legendLab1) +
    annotation_custom(legendLab2) +
    theme_bw() +
    theme(
          axis.ticks = element_line(size = 1.0, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size = 22, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 30, colour = "black"),
          legend.position = "none",
          #legend.text = element_text(size = 10),
          #legend.background = element_rect(fill = "transparent"),
          #legend.key = element_rect(colour = "transparent",
          #                          fill = "transparent"),
          #legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 3.5, colour = "black"),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote("Random loci" ~
                   "(" * italic("n") ~ "=" ~
                   .(prettyNum(summaryDFranLoc2$n[1],
                               big.mark = ",", trim = T)) *
                   ")"))
  ggObjGA <- grid.arrange(ggObj1, ggObj2, ggObj3, nrow = 1, ncol = 3)
  ggsave(paste0(plotDir,
                paste(libNames, collapse = "_"),
                "_around_", featureName, "_", superfamNames[y], ".pdf"),
         plot = ggObjGA,
         height = 6.5, width = 21)
}
