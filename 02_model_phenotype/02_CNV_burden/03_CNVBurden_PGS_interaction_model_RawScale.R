# Model raw phenotypes according to CNV burden, PGS, and their interaction for 43 assessed phenotypes

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

# Phenotype (raw scale; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing raw phenotype values 
pheno <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_WB_raw.txt"))
pheno <- pheno[pheno$IID %in% test_samples$IID, ]

# Polygenic scores (PGSs; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous.txt.gz"))
pgs <- pgs[pgs$IID %in% test_samples$IID, ]

# Covariates (filtered for testing samples IID) )
# File with sample identifier (IID) as first column, then one column per covariate, including age, age^2, sex, genotyping batch, and PC1-40
cov <- as.data.frame(fread("CNV_PGS/data/covariates.txt"))
cov <- cov[cov$IID %in% test_samples$IID, ]

# CNV burden data (filtered for testing samples IID)
# CNV burdens for the UKBB (as described in Auwerx et al., 2022)
cnv_burden <- as.data.frame(fread("CNV_PGS/data/CNV_burden.txt.gz", select = c(1,10), col.names = c("IID", "CNV_GENES")))
dup_burden <- as.data.frame(fread("CNV_PGS/data/DUP_burden.txt.gz", select = c(1,10), col.names = c("IID", "DUP_GENES")))
del_burden <- as.data.frame(fread("CNV_PGS/data/DEL_burden.txt.gz", select = c(1,10), col.names = c("IID", "DEL_GENES")))
cnv <- left_join(cnv_burden, dup_burden, by = "IID")
cnv <- left_join(cnv, del_burden, by = "IID")
cnv <- cnv[cnv$IID %in% test_samples$IID, ]
rm(cnv_burden, dup_burden, del_burden)


########################################################
# STEP 2: Regression analysis
########################################################

# Create a dataframe to store results
df <- cnv_signals
rm(cnv_signals)

# Loop over phenotypes
for (i in 1:nrow(df)) {
  
  # Define the phenotype
  p <- df[i, "PHENO"]
  
  # Merge phenotype, PGS, CNV burden, and covariates
  df_temp <- left_join(pheno[, names(pheno) %in% c("IID", p)], cnv, by = "IID")
  colnames(df_temp) [2] <- "PHENO"
  df_temp <- left_join(df_temp, pgs[, names(pgs) %in% c("IID", p)], by = "IID")
  colnames(df_temp) [6] <- "PGS"
  df_temp <- na.omit(left_join(df_temp, cov, by = "IID"))

  # Fit linear regression - CNV burden
  fit_CNV <- lm(PHENO ~ . + CNV_GENES:PGS, data = df_temp[,-c(1,4,5)])
  
  # Fill in the result table
  df[i, "EFFECT_CNV_GENES"] <- summary(fit_CNV)$coefficients[2,1]
  df[i, "SE_CNV_GENES"] <- summary(fit_CNV)$coefficients[2,2]
  df[i, "P_CNV_GENES"] <- summary(fit_CNV)$coefficients[2,4]
  
  df[i, "EFFECT_cnv_PGS"] <- summary(fit_CNV)$coefficients[3,1]
  df[i, "SE_cnv_PGS"] <- summary(fit_CNV)$coefficients[3,2]
  df[i, "P_cnv_PGS"] <- summary(fit_CNV)$coefficients[3,4]
  
  df[i, "EFFECT_CNV_GENESxPGS"] <- summary(fit_CNV)$coefficients[152,1]
  df[i, "SE_CNV_GENESxPGS"] <- summary(fit_CNV)$coefficients[152,2]
  df[i, "P_CNV_GENESxPGS"] <- summary(fit_CNV)$coefficients[152,4]
  
  df[i, "adjR2_CNV_GENES"] <- summary(fit_CNV)$adj.r.squared
   
  # Fit linear regression - DEL burden
  fit_DEL <- lm(PHENO ~ . + DEL_GENES:PGS, data = df_temp[,-c(1,3,4)])
  
  # Fill in the result table
  df[i, "EFFECT_DEL_GENES"] <- summary(fit_DEL)$coefficients[2,1]
  df[i, "SE_DEL_GENES"] <- summary(fit_DEL)$coefficients[2,2]
  df[i, "P_DEL_GENES"] <- summary(fit_DEL)$coefficients[2,4]
  
  df[i, "EFFECT_del_PGS"] <- summary(fit_DEL)$coefficients[3,1]
  df[i, "SE_del_PGS"] <- summary(fit_DEL)$coefficients[3,2]
  df[i, "P_del_PGS"] <- summary(fit_DEL)$coefficients[3,4]
  
  df[i, "EFFECT_DEL_GENESxPGS"] <- summary(fit_DEL)$coefficients[152,1]
  df[i, "SE_DEL_GENESxPGS"] <- summary(fit_DEL)$coefficients[152,2]
  df[i, "P_DEL_GENESxPGS"] <- summary(fit_DEL)$coefficients[152,4]
  
  df[i, "adjR2_DEL_GENES"] <- summary(fit_DEL)$adj.r.squared
  
  # Fit linear regression - DUP burden
  fit_DUP <- lm(PHENO ~ . + DUP_GENES:PGS, data = df_temp[,-c(1,3,5)])
  
  # Fill in the result table
  df[i, "EFFECT_DUP_GENES"] <- summary(fit_DUP)$coefficients[2,1]
  df[i, "SE_DUP_GENES"] <- summary(fit_DUP)$coefficients[2,2]
  df[i, "P_DUP_GENES"] <- summary(fit_DUP)$coefficients[2,4]
  
  df[i, "EFFECT_dup_PGS"] <- summary(fit_DUP)$coefficients[3,1]
  df[i, "SE_dup_PGS"] <- summary(fit_DUP)$coefficients[3,2]
  df[i, "P_dup_PGS"] <- summary(fit_DUP)$coefficients[3,4]
  
  df[i, "EFFECT_DUP_GENESxPGS"] <- summary(fit_DUP)$coefficients[152,1]
  df[i, "SE_DUP_GENESxPGS"] <- summary(fit_DUP)$coefficients[152,2]
  df[i, "P_DUP_GENESxPGS"] <- summary(fit_DUP)$coefficients[152,4]
  
  df[i, "adjR2_DUP_GENES"] <- summary(fit_DUP)$adj.r.squared
  
  }
rm(i, p, df_temp, fit_CNV, fit_DUP, fit_DEL)


########################################################
# STEP 3: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_phenotype/CNVBurden_PGS_interaction_model_RawScale.txt", col.names = T, row.names = F, quote = F, sep = "\t")
