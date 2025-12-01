# Model adjusted phenotypes according to CNV status and PGS for 119 CNV-trait pairs

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
# File with a single column (IID) containing the sample identifier for all samples in the test set. 
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# CNV signals
# File with the following columns: PHENO, CHR, CNVR_START, CNVR_STOP, TOP_MODEL, CB (cytogenic band)
cnv_signals <- as.data.frame(fread("CNV_PGS/data/cnv_signals.txt"))

# Phenotypic variability explained by the PGS
# File with the following columns:  PHENO, cor_2 (squared correlation between phenotype and PGS)
cor2 <- as.data.frame(fread("CNV_PGS/data/phenotype_explained_by_PGS.txt"))

# Phenotype (covariate-corrected + INT; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing covariate-adjusted, inverse-normal transformed phenotype values 
pheno <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.txt"))
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
# STEP 2: Regression analysis
########################################################

# Create a dataframe to store results
df <- right_join(cor2, cnv_signals, by = "PHENO")
rm(cor2)

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

  # Identify high confidence deletion and duplication carriers
  del_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  dup_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  
  # Merge Pheno and PGS data
  df_temp <- left_join(pheno[, names(pheno) %in% c("IID", p)], pgs[, names(pgs) %in% c("IID", p)], by = "IID")
  colnames(df_temp) <- c("IID", "PHENO", "PGS")

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

  # Fit linear regression
  fit <- lm(PHENO ~ CNV + PGS, data = df_temp)
  
  # Fill in the result table with sample number and regression coefficients
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
  
  df[i, "EFFECT_CNV"] <- summary(fit)$coefficients[2,1]
  df[i, "SE_CNV"] <- summary(fit)$coefficients[2,2]
  df[i, "P_CNV"] <- summary(fit)$coefficients[2,4]
  
  df[i, "EFFECT_PGS"] <- summary(fit)$coefficients[3,1]
  df[i, "SE_PGS"] <- summary(fit)$coefficients[3,2]
  df[i, "P_PGS"] <- summary(fit)$coefficients[3,4]
  
  df[i, "adjR2"] <- summary(fit)$adj.r.squared
  
}
rm(i, p, chr, pos, model, del_carriers, dup_carriers, df_temp)


########################################################
# STEP 3: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_phenotype/CNV_PGS_model.txt", col.names = T, row.names = F, quote = F, sep = "\t")
