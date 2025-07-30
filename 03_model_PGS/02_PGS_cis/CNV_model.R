# Model local PGS_cis according to CNV status for 119 CNV-trait pairs

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

# Covariates (raw; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then age, age2^2, sex, array (i.e., genotyping array), batch (i.e., genotyping batch), principal components 1-40
cov <- as.data.frame(fread("CNV_PGS/data/covariates.txt"))
cov <- cov[cov$IID %in% test_samples$IID, ]

# PGS data (50kb in) for all individuals (local PGS_cis ± 50kb around CNV region; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous.in.50kb.txt.gz"))
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

# Create a dataframe to store regression results
df <- right_join(cor2, cnv_signals, by = "PHENO")
rm(cor2)

# Loop over CNV-trait pairs
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
  
  # Subset PGS data
  df_temp <- pgs[, names(pgs) %in% c("IID", paste0(p, "_chr", chr, "_", start-50000, "_", end+50000))]
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
  
  # Add covariates
  df_temp <- left_join(df_temp, cov[, -c(2:6, 27:46)], by = "IID")
  
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

  # Fill in the table with correlation coefficients
  df[i, "cor_PEARSON"] <- cor.test(df_temp$PGS, df_temp$CNV, method = "pearson")$estimate
  df[i, "P_PEARSON"] <- cor.test(df_temp$PGS, df_temp$CNV, method = "pearson")$p.val
  
  # Fit linear regression
  fit <- lm(PGS ~ CNV, data = df_temp[, -c(1)])
  
  # Fill table with results
  df[i, "EFFECT_CNV_PGS_cis"] <- summary(fit)$coefficients[2,1]
  df[i, "SE_CNV_PGS_cis"] <- summary(fit)$coefficients[2,2]
  df[i, "P_CNV_PGS_cis"] <- summary(fit)$coefficients[2,4]
  df[i, "adjR2_CNV_PGS_cis"] <- summary(fit)$adj.r.squared
  
}
rm(i, p, chr, pos, model, del_carriers, dup_carriers, df_temp)


########################################################
# STEP 3: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_PGS/CNV_model_PGS_cis.txt", col.names = T, row.names = F, quote = F, sep = "\t")
