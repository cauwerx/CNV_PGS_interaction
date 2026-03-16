# Model genome-wide PGS_gw according to CNV status for 119 CNV-trait pairs

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)


########################################################
# STEP 1: Load data
########################################################

# Testing samples
# File with a single column (IID) containing the sample identifier for all samples in the test set
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# CNV signals
# File with the following columns: PHENO, CHR, CNVR_START, CNVR_STOP, TOP_MODEL, CB (cytogenic band)
cnv_signals <- as.data.frame(fread("CNV_PGS/data/cnv_signals.txt"))

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
# STEP 2: Regression analysis
########################################################

# Create a dataframe to store results
df <- cnv_signals
rm(cnv_signals)

# Loop over CNV-trait pairs
for (i in 1:nrow(df)) {
  
  # Define the signal
  p <- df[i, "PHENO"]
  chr <- df[i, "CHR"]
  pos <-  df[i, "TOP_POS"]
  model <- df[i, "TOP_MODEL"]; if(model == "M-DUP") {model <- "M"}
  print(paste0("Analyzing ", p," for chr", chr, ":", pos, " CNVs (", model, ")"))
  
  # Identify samples carrying relevant CNVs (overlapping the lead probe)
  df_cnvs <- cnvs[which(cnvs$Chromosome == chr & cnvs$Start_Position_bp <= pos & cnvs$End_Position_bp >= pos), ]

  # Identify high-confidence deletion and duplication carriers
  del_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  dup_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  
  # Subset PGS data
  df_temp <- pgs[, names(pgs) %in% c("IID", p)]
  colnames(df_temp) <- c("IID", "PGS")

  # Add CNV carrier status for the relevant model
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

  # Fill in the result table with sample counts
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
  
  # Fit linear regression
  fit <- lm(PGS ~ CNV, data = df_temp[, -c(1)])
  
  # Fill table with results
  df[i, "EFFECT_CNV_PGS_GW"] <- summary(fit)$coefficients[2,1]
  df[i, "SE_CNV_PGS_GW"] <- summary(fit)$coefficients[2,2]
  df[i, "P_CNV_PGS_GW"] <- summary(fit)$coefficients[2,4]
  df[i, "adjR2_CNV_PGS_GW"] <- summary(fit)$adj.r.squared

}
rm(i, p, chr, pos, model, del_carriers, dup_carriers, df_temp)


########################################################
# STEP 3: Compare results to CNV effect on phenotype
########################################################

# Load CNV effects on phenotype; Output of 02_model_phenotype/01_CNV/01_CNV_model.R
cnv_effect_pheno <- as.data.frame(fread("CNV_PGS/data/model_phenotype/CNV_model.txt"))
colnames(cnv_effect_pheno)[c(13:15)] <- c("EFFECT_CNV_PHENO", "SE_CNV_PHENO", "P_CNV_PHENO")

# Merge the CNV effect on PGS to the CNV effect on the phenotype
cnv_effect <- left_join(df, cnv_effect_pheno)

# Sign concordance (i.e., if CNV increases a trait, the CNV will tend to increase the PGS for that trait & vice versa)

# Among all CNV-trait pairs
nrow(cnv_effect[which(cnv_effect$EFFECT_CNV_PGS_GW > 0 & cnv_effect$EFFECT_CNV_PHENO > 0), ]) + nrow(cnv_effect[which(cnv_effect$EFFECT_CNV_PGS_GW < 0 & cnv_effect$EFFECT_CNV_PHENO < 0), ])

# Among nominally significant 
nrow(cnv_effect[which(cnv_effect$EFFECT_CNV_PGS_GW > 0 & cnv_effect$EFFECT_CNV_PHENO > 0 & cnv_effect$P_CNV_PGS_GW < 0.05), ]) + nrow(cnv_effect[which(cnv_effect$EFFECT_CNV_PGS_GW < 0 & cnv_effect$EFFECT_CNV_PHENO < 0 & cnv_effect$P_CNV_PGS_GW < 0.05), ])


########################################################
# STEP 4: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_PGS/CNV_model_PGS_GW.txt", col.names = T, row.names = F, quote = F, sep = "\t")
