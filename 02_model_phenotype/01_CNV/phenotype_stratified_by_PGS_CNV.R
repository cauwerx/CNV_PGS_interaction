# Raw phenotypic measurements stratified by sample CNV carrier status and PGS quintile for 119 CNV-trait pairs

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(MetBrewer)
library(ggplot2)
library(showtext)
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")


########################################################
# STEP 1: Load data
########################################################

# Testing samples; 
# File with a single column (IID) containing the sample identifier for all samples in the test set. 
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# CNV signals
# File with the following columns: PHENO, CHR, CNVR_START, CNVR_STOP, TOP_MODEL, CB (cytogenic band), Y_AXIS (contains the label for the y-axis, as desired)
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
# STEP 2: Stratify phenotype by CNV carrier status and PGS
########################################################

# Create a list to store the results
plot_list <- list()

# Loop over CNV-trait pairs
for (i in 1:nrow(cnv_signals)) {
   
    # Define the signal
    p <- cnv_signals[i, "PHENO"]
    chr <- cnv_signals[i, "CHR"]
    pos <-  cnv_signals[i, "TOP_POS"]
    model <- cnv_signals[i, "TOP_MODEL"]
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

    # DELETION: Merge all data into a single dataframe; split into quintiles (providing there are > 15 DEL carriers)
    if(nrow(df_del) > 0) {
      df_del <- left_join(df_del[, names(df_del) %in% c("IID", "TYPE")], pheno[, names(pheno) %in% c("IID", p)], by = "IID")
      colnames(df_del)[3] <- "PHENO"
    }
    if(nrow(df_del) >= 15) {
       df_del <- left_join(df_del, pgs[, names(pgs) %in% c("IID", p)], by = "IID")
       colnames(df_del)[4] <- "PGS"
       df_del <- na.omit(df_del)
       df_del$TYPE <- paste0("Deletion\n(N = ", nrow(df_del), ")")
       quintile_breaks_del <- quantile(df_del$PGS, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
       df_del$Quintile <- cut(df_del$PGS, breaks = quintile_breaks_del, labels = FALSE, include.lowest = TRUE)
       rm(quintile_breaks_del)
       df_del$GROUP <- factor(df_del$Quintile, levels = c(1,2,3,4,5), labels = c("1st","2nd","3rd","4th","5th"))
    }
     
    # DUPLICATION: Merge all data into a single dataframe; split into quintiles (providing there are > 15 DUP carriers)
    if(nrow(df_dup) > 0) {
      df_dup <- left_join(df_dup[, names(df_dup) %in% c("IID", "TYPE")], pheno[, names(pheno) %in% c("IID", p)], by = "IID")
      colnames(df_dup)[3] <- "PHENO"
    }
    if(nrow(df_dup) >= 15) {
       df_dup <- left_join(df_dup, pgs[, names(pgs) %in% c("IID", p)], by = "IID")
       colnames(df_dup)[4] <- "PGS"
       df_dup <- na.omit(df_dup)
       df_dup$TYPE <- paste0("Duplication\n(N = ", nrow(df_dup), ")")
       quintile_breaks_dup <- quantile(df_dup$PGS, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
       df_dup$Quintile <- cut(df_dup$PGS, breaks = quintile_breaks_dup, labels = FALSE, include.lowest = TRUE)
       rm(quintile_breaks_dup)
       df_dup$GROUP <- factor(df_dup$Quintile, levels = c(1,2,3,4,5), labels = c("1st","2nd","3rd","4th","5th"))
    }
    
    # COPY-NEUTRAL: Merge all data into a single dataframe; split into quintiles
    df_cn <- left_join(df_cn, pheno[, names(pheno) %in% c("IID", p)], by = "IID")
    colnames(df_cn)[3] <- "PHENO"
    df_cn <- left_join(df_cn, pgs[, names(pgs) %in% c("IID", p)], by = "IID")
    colnames(df_cn)[4] <- "PGS"
    df_cn <- na.omit(df_cn)
    df_cn$TYPE <- paste0("Copy-neutral\n(N = ", nrow(df_cn), ")")
    quintile_breaks_cn <- quantile(df_cn$PGS, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
    df_cn$Quintile <- cut(df_cn$PGS, breaks = quintile_breaks_cn, labels = FALSE, include.lowest = TRUE)
    rm(quintile_breaks_cn)
    df_cn$GROUP <- factor(df_cn$Quintile, levels = c(1,2,3,4,5), labels = c("1st","2nd","3rd","4th","5th"))
    
    # Merge the datasets (among DEL, DUP, and copy-neutral) with > 15 samples
    if(nrow(df_del) >= 15 & nrow(df_dup) >= 15) {df <- rbind(df_del, df_dup)}
    if(nrow(df_del) >= 15 & nrow(df_dup) < 15) {df <- df_del}
    if(nrow(df_del) < 15 & nrow(df_dup) >= 15) {df <- df_dup}
    if(nrow(df_del) < 15 & nrow(df_dup) < 15) {next}
    df <- rbind(df, df_cn)
  
    # Factorize
    df$TYPE <- factor(df$TYPE, levels =  unique(df$TYPE))
    if(length(unique(df$TYPE)) == 3) {
      df$TYPE <- factor(df$TYPE, levels = unique(df$TYPE)[c(1,3,2)])
    } else if (any(grepl("^Del", unique(df$TYPE)))){
      df$TYPE <- factor(df$TYPE, levels = unique(df$TYPE)[c(1,2)])
    } else if (any(grepl("^Dup", unique(df$TYPE)))) {
      df$TYPE <- factor(df$TYPE, levels = unique(df$TYPE)[c(2,1)])
    }
    
    # Add Header
    df$HEADER <- factor(cnv_signals[i, "CB"])
   
    # Calculate y-limit limits for plotting
    y_lim_range <- vector()
    for (t in unique(df$TYPE)) {
      for (q in unique(df$Quintile)) {
        y_lim_range <- c(y_lim_range, boxplot.stats(df$PHENO[df$TYPE %in% t & df$Quintile %in% q])$stats[c(1,5)])
      }
    }; rm(t)
    y_lim_range <- range(y_lim_range)

    # Calculate median phenotype
    median_pheno <- median(pheno[!pheno$IID %in% c(df_cnvs$IID), p], na.rm = T)
    rm(df_cnvs, df_del, df_dup, df_cn)
                 
    # Plot
    plot_list[[i]] <- ggplot(df) +
                        # Facet
                        facet_wrap(.~ HEADER) +
                        # Add median phenotype value
                        geom_hline(yintercept = median_pheno, color = "brown4", linewidth = 1.3, alpha = 0.7) +
                        # Boxplot
                        geom_boxplot(aes(x = TYPE, y = PHENO, fill = GROUP, alpha = GROUP), outlier.shape = NA, position = position_dodge(width = 0.85), width = 0.6, fill = "#2B4050") + 
                        scale_fill_manual(paste0(p, "\nPGS Quintile"), values = c("#2B4050", "#2B4050", "#2B4050",  "#2B4050", "#2B4050")) +
                        scale_alpha_manual(paste0(p, "\nPGS Quintile"), values = c(0.15, 0.35, 0.55, 0.75, 0.95)) +
                        # Layout
                        xlab("") +  ylab(cnv_signals[i, "Y_AXIS"]) +
                         coord_cartesian(ylim = y_lim_range + c(-0.2, 0.2) * diff(y_lim_range) / 2) +
                         theme_bw() +
                         theme(text = element_text(family = "Arial"),
                               axis.text = element_text(size = 11),
                               axis.title.y = element_text(size = 11, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                               legend.title = element_text(size = 11, face = "bold"),
                               strip.text = element_text(size = 12, face = "bold", color = "#2B4050"),
                               strip.background = element_rect(fill = "#F4F4F4", color = "#2B4050"))
    
    # Plot individual data points only if there are less than 1000 DEL carreiers
    if(nrow(df[grep("^Del", df$TYPE), ]) < 1000) {
      plot_list[[i]] <- plot_list[[i]] + geom_point(df[grep("^Del", df$TYPE), ], mapping = aes(x = TYPE, y = PHENO, fill = GROUP, alpha = GROUP), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85), alpha = 0.25, color = "#2B4050", size = 1)
    }
    # Plot individual data points only if there are less than 1000 DUP carreiers
    if(nrow(df[grep("^Dup", df$TYPE), ]) < 1000) {
      plot_list[[i]] <- plot_list[[i]] + geom_point(df[grep("^Dup", df$TYPE), ], mapping = aes(x = TYPE, y = PHENO, fill = GROUP, alpha = GROUP), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85), alpha = 0.25, color = "#2B4050", size = 1)
    }
}
rm(i, p, chr, pos, model, y_lim_range, median_pheno)


########################################################
# STEP 3: Save plots
########################################################

# Save all plots in a single PDF
pdf("CNV_PGS/plots/phenotype_stratified_PGS_CNV.pdf", width = 6, height = 2.5)
showtext_auto() 
for (p in plot_list) {
  print(p)
}
dev.off()
showtext_auto(FALSE) 