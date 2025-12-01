# PGSs stratified by sample CNV carrier status for 119 CNV-trait pairs

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)
library(MetBrewer)
library(ggplot2)
library(ggsignif)
library(lazyWeave)
library(ggpubr)
library(showtext)
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")


########################################################
# STEP 1: Load data
########################################################

# Testing samples; 
# File with a single column (IID) containing the sample identifier for all samples in the test set. 
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# CNV signals
# File with the following columns: PHENO, CHR, CNVR_START, CNVR_STOP, TOP_MODEL, CB (cytogenic band)
cnv_signals <- as.data.frame(fread("CNV_PGS/data/cnv_signals.txt"))

# Phenotype (raw; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing raw phenotype values 
pheno <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_WB_raw.txt"))
pheno <- pheno[pheno$IID %in% test_samples$IID, ]

# Polygenic scores (PGSs; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous.txt.gz"))
pgs <- pgs[pgs$IID %in% test_samples$IID, ]

# CNV calls (filtered for testing samples IID)
# CNV calls for the UKBB (as described in Auwerx et al., 2022)
cnvs <- fread("CNV_PGS/data/ukb_cnv_global.gz", header = T)
cnvs <- separate(cnvs, Sample_Name, c("IID", "batch"), sep = "_", remove = FALSE, extra = "merge")
cnvs$IID <- as.integer(cnvs$IID)
cnvs <- as.data.frame(cnvs[cnvs$IID %in% test_samples$IID, ])


########################################################
# STEP 2: Stratify PGS by CNV carrier status
########################################################

# Create a list to store the results
plot_list <- list()

# Loop over CNV-trait pairs
for (i in 1:nrow(cnv_signals)) {
  
  # Define the signal
  p <- cnv_signals[i, "PHENO"]
  chr <- cnv_signals[i, "CHR"]
  pos <-  cnv_signals[i, "TOP_POS"]
  if(cnv_signals[i, "TOP_MODEL"] == "M-DUP") {model <- "DUP"} else {model <- cnv_signals[i, "TOP_MODEL"]}
  if(cnv_signals[i, "P_CNV_PGS_GW"] > 0.05) {
    pval_label <-  paste0("p[", model, "] == ns")
  } else {
      pval_label <- paste0("p[", model, "] == ", pvalString(cnv_signals[i, "P_CNV_PGS_GW"], format = "exact", digits = 2))
  }
  print(paste0("Analyzing ", p," among carriers of chr", chr, ":", pos, " CNVs (", model, ")"))
  
  # Identify samples carrying relevant CNVs (overlapping the lead signal)
  df_cnvs <- cnvs[which(cnvs$Chromosome == chr & cnvs$Start_Position_bp <= pos & cnvs$End_Position_bp >= pos), ]
  
  # Split CNV carriers by CNV type
  df_del <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), ]
  if(nrow(df_del) > 0) {df_del$TYPE <- "Deletion"}
  df_dup <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), ]
  if(nrow(df_dup) > 0) {df_dup$TYPE <- "Duplication"}
  
  # Get the copy-neutral individuals (no CNV, not even with a low quality score)
  df_cn <- pheno[, "IID", drop = F]
  df_cn <- df_cn[!df_cn$IID %in% df_cnvs$IID, , drop = F]
  df_cn$TYPE <- "Copy-Neutral"
  
  # DELETION: Merge all data into a single dataframe
  if(nrow(df_del) > 0) {
    df_del <- na.omit(left_join(df_del[, names(df_del) %in% c("IID", "TYPE")], pgs[, names(pgs) %in% c("IID", p)], by = "IID"))
    colnames(df_del)[3] <- "PGS"
    df_del$TYPE <- paste0("Deletion\n(N = ", nrow(df_del), ")")
  }

  # DUPLICATION: Merge all data into a single dataframe
  if(nrow(df_dup) > 0) {
    df_dup <- na.omit(left_join(df_dup[, names(df_dup) %in% c("IID", "TYPE")], pgs[, names(pgs) %in% c("IID", p)], by = "IID"))
    colnames(df_dup)[3] <- "PGS"
    df_dup$TYPE <- paste0("Duplication\n(N = ", nrow(df_dup), ")")
  }
  
  # COPY-NEUTRAL: Merge all data into a single dataframe
  df_cn <- na.omit(left_join(df_cn, pgs[, names(pgs) %in% c("IID", p)], by = "IID"))
  colnames(df_cn)[3] <- "PGS"
  df_cn$TYPE <- paste0("Copy-neutral\n(N = ", nrow(df_cn), ")")

  # Merge the datasets (among DEL, DUP, and copy-neutral) with > 0 samples
  if(nrow(df_del) > 0 & nrow(df_dup) > 0) {df <- rbind(df_del, df_dup)}
  if(nrow(df_del) > 0 & nrow(df_dup) == 0) {df <- df_del}
  if(nrow(df_del) == 0 & nrow(df_dup) > 0) {df <- df_dup}
  df <- rbind(df, df_cn)
  
  # Factorize
  if(length(unique(df$TYPE)) == 3) {
    df$TYPE <- factor(df$TYPE, levels = unique(df$TYPE)[c(1,3,2)])
  } else if (any(grepl("^Del", unique(df$TYPE)))){
    df$TYPE <- factor(df$TYPE, levels = c(unique(df$TYPE)[c(1,2)], "Duplication\n(N = 0)"))
  } else if (any(grepl("^Dup", unique(df$TYPE)))) {
    df$TYPE <- factor(df$TYPE, levels = c("Deletion\n(N = 0)", unique(df$TYPE)[c(2,1)]))
 
  }
  # Add Header
  df$HEADER <- factor(paste0(cnv_signals[i, "CB"]))
  
  # Calculate y-limit limits for plotting
  y_lim_range <- vector()
  for (t in unique(df$TYPE)) {
      y_lim_range <- c(y_lim_range, boxplot.stats(df$PGS[df$TYPE %in% t])$stats[c(1,5)])
  }; rm(t)
  y_lim_range <- range(y_lim_range)
  
  # Significance table
  if(model == "M") {
    sig_df <- data.frame(start = factor(levels(df$TYPE)[grep("Del", levels(df$TYPE))]), 
                         end = factor(levels(df$TYPE)[grep("Dup", levels(df$TYPE))]), 
                         y = y_lim_range[2] + 0.2 * diff(y_lim_range) / 2,
                         label = pval_label)
  } else if (model == "DEL") {
    sig_df <- data.frame(start = factor(levels(df$TYPE)[grep("Del", levels(df$TYPE))]), 
                         end = factor(levels(df$TYPE)[grep("Copy", levels(df$TYPE))]), 
                         y = y_lim_range[2] + 0.2 * diff(y_lim_range) / 2,
                         label = pval_label)
  }  else if (model == "DUP") {
    sig_df <- data.frame(start = factor(levels(df$TYPE)[grep("Copy", levels(df$TYPE))]), 
                         end = factor(levels(df$TYPE)[grep("Dup", levels(df$TYPE))]), 
                         y = y_lim_range[2] + 0.2 * diff(y_lim_range) / 2,
                         label = pval_label)
  
  # Plot
  plot_list[[i]] <- ggplot(df) +
                    # Facet
                    facet_wrap(.~ HEADER) +
                    # Boxplot
                    geom_boxplot(aes(x = TYPE, y = PGS, fill = TYPE), outlier.shape = NA, position = position_dodge(width = 0.85), width = 0.6) +
                    scale_fill_manual("CNV status", values = c("indianred2", "#2B4050", "cornflowerblue"), drop = FALSE, guide = "none") +
                    # P-value
                    geom_signif(data = sig_df, aes(xmin = start, xmax = end, annotations = label, y_position = y), textsize = 3.5, vjust = -0.4, manual = TRUE, tip_length = 0, parse = T) +
                    # Layout
                    ylab(paste0(label, " PGS")) +
                    scale_x_discrete(labels = levels(df$TYPE), drop = F) +
                    coord_cartesian(ylim = y_lim_range + c(-0.15, 0.5) * diff(y_lim_range) / 2) +
                    theme_bw() +
                    theme(text = element_text(family = "Arial"),
                          axis.text.x = element_text(size = 9),
                          axis.text.y = element_text(size = 11),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size = 11, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                          strip.text = element_text(size = 12, face = "bold", color = "#2B4050"),
                          strip.background = element_rect(fill = "#F4F4F4", color = "#2B4050"))
  
    # Plot individual data points only if there are less than 1000 DEL carreiers
    if(nrow(df[grep("^Del", df$TYPE), ]) < 1000) {
      plot_list[[i]] <- plot_list[[i]] + geom_point(df[grep("^Del", df$TYPE), ], mapping = aes(x = TYPE, y = PGS, fill = TYPE), position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.85), alpha = 0.25, color = "#2B4050", size = 1)
    }
  
    # Plot individual data points only if there are less than 1000 DUP carreiers
    if(nrow(df[grep("^Dup", df$TYPE), ]) < 1000) {
      plot_list[[i]] <- plot_list[[i]] + geom_point(df[grep("^Dup", df$TYPE), ], mapping = aes(x = TYPE, y = PGS, fill = TYPE), position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.85), alpha = 0.25, color = "#2B4050", size = 1)
    }                
}
rm(i, p, chr, pos, model, pval_label, df_cn, df_cnvs, df_del, df_dup, y_lim_range, sig_df)


########################################################
# STEP 3: Save plots
########################################################

# Save all plots in a single PDF
pdf("CNV_PGS/plots/PGS_stratified_CNV.pdf", width = 3.2, height = 2.45)
showtext_auto() 
for (p in plot_list) {
  print(p)
}
dev.off()
showtext_auto(FALSE)
