# Correlation between PGS_cis [CHR] and CNV carrier status for 119 CNV-trait pairs

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)


########################################################
# STEP 1: Load data
########################################################

# Testing samples; 
# File with a single column (IID) containing the sample identifier for all samples in the test set 
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# CNV signals
# File with the following columns: PHENO, CHR, CNVR_START, CNVR_STOP, TOP_MODEL, CB (cytogenic band)
cnv_signals <- as.data.frame(fread("CNV_PGS/data/cnv_signals.txt"))

# PGS data for all individuals (PGS_cis ± chromosome around CNV region; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous_in_CHR.txt.gz"))
pgs <- pgs[pgs$IID %in% test_samples$IID, ]

# CNV calls (filtered for testing samples IID)
# CNV calls for the UKBB (as described in Auwerx et al., 2022)
cnvs <- fread("CNV_PGS/data/ukb_cnv_global.gz", header = T)
cnvs <- separate(cnvs, Sample_Name, c("IID", "batch"), sep = "_", remove = FALSE, extra = "merge")
cnvs$IID <- as.integer(cnvs$IID)
cnvs <- as.data.frame(cnvs[cnvs$IID %in% test_samples$IID, ])


########################################################
# STEP 2: Correlation analysis
########################################################

# Create a dataframe to store correlation results
df <- cnv_signals
rm(cnv_signals)

# Loop over phenotypes
for (i in 1:nrow(df)) {
  
  # Define the signal
  p <- df[i, "PHENO"]
  chr <- df[i, "CHR"]
  pos <-  df[i, "TOP_POS"]
  start <- df[i, "CNVR_START"]
  end <- df[i, "CNVR_STOP"]
  model <- df[i, "TOP_MODEL"]; if(model == "M-DUP") {model <- "M"}
  print(paste0("Analyzing ", p," for chr", chr, ":", pos, " CNVs (", model, ")"))
  
  # Identify samples carrying relevant CNVs
  df_cnvs <- cnvs[which(cnvs$Chromosome == chr & cnvs$Start_Position_bp <= pos & cnvs$End_Position_bp >= pos), ]
  
  # Identify high confidence deletion and duplication carriers
  del_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  dup_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  
  # Merge PGS data & add CNV carrier status for the relevant model
  df_temp <- pgs[, names(pgs) %in% c("IID", paste0(p, "_chr", chr))]
  colnames(df_temp) <- c("IID", "PGS")
  df_temp$CNV <- 0
  df_temp[df_temp$IID %in% df_cnvs$IID, "CNV"] <- NA
  if(model == "DEL") {
    df_temp[df_temp$IID %in% del_carriers, "CNV"] <- 1
  }
  if(model == "DUP") {
    df_temp[df_temp$IID %in% dup_carriers, "CNV"] <- 1
  }
  if(model == "M") {
    df_temp[df_temp$IID %in% del_carriers, "CNV"] <- -1
    df_temp[df_temp$IID %in% dup_carriers, "CNV"] <- 1
  }
  df_temp <- na.omit(df_temp)
  
  # Get CNV carrier counts
  if(model == "DEL") {
    df[i, "N_DEL"] <- table(df_temp$CNV, useNA = "always")[2]
    df[i, "N_DUP"] <- NA
    df[i, "N_COPY_NEUTRAL"] <- table(df_temp$CNV, useNA = "always")[1]
  }
  if(model == "DUP") {
    df[i, "N_DEL"] <- NA
    df[i, "N_DUP"] <- table(df_temp$CNV, useNA = "always")[2]
    df[i, "N_COPY_NEUTRAL"] <- table(df_temp$CNV, useNA = "always")[1]
  }
  if(model == "M") {
    df[i, "N_DEL"] <- table(df_temp$CNV, useNA = "always")[1]
    df[i, "N_DUP"] <- table(df_temp$CNV, useNA = "always")[3]
    df[i, "N_COPY_NEUTRAL"] <- table(df_temp$CNV, useNA = "always")[2]
  }
  
  # Get correlation coefficients
  df[i, "cor_PEARSON"] <- cor.test(df_temp$PGS, df_temp$CNV, method = "pearson")$estimate
  df[i, "P_PEARSON"] <- cor.test(df_temp$PGS, df_temp$CNV, method = "pearson")$p.val
  
  # Get correlation 95%-CI by Fisher's z-transformation
  z_r <- 0.5 * log((1 + df[i, "cor_PEARSON"]) / (1 - df[i, "cor_PEARSON"]))
  se_z <- 1 / sqrt((df[i, "N_DEL"] + df[i, "N_DUP"] + df[i, "N_COPY_NEUTRAL"]) - 3)                   
  z_L95 <- z_r - (1.96 * se_z)
  z_U95 <- z_r + (1.96 * se_z)
  df[i, "cor_PEARSON_L95"] <- (exp(2 * z_L95) - 1) / (exp(2 * z_L95) + 1)
  df[i, "cor_PEARSON_U95"] <- (exp(2 * z_U95) - 1) / (exp(2 * z_U95) + 1)
  
}
rm(i, p, chr, pos, model, del_carriers, dup_carriers, df_temp, z_r, se_z, z_L95, z_U95)


########################################################
# STEP 3: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_PGS/CNV_cor_PGS_cis_CHR.txt", col.names = T, row.names = F, quote = F, sep = "\t")
