# Model adjusted phenotypes according to CNV burden for 43 assessed phenotypes

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

# Phenotypic variability explained by the PGS
# File with the following columns:  PHENO, cor_2 (squared correlation between phenotype and PGS)
cor2 <- as.data.frame(fread("CNV_PGS/data/phenotype_explained_by_PGS.txt"))

# Phenotype (covariate-corrected + INT; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing covariate-adjusted, inverse-normal transformed phenotype values 
pheno <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.txt"))
pheno <- pheno[pheno$IID %in% test_samples$IID, ]

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
# STEP 2:Regression analysis
########################################################

# Create a dataframe to store regression results
df <- cor2[, c(1,5)]
rm(cor2)

# Loop over phenotypes
for (i in 1:nrow(df)) {
  
  # Define the signal
  p <- df[i, "PHENO"]

  # Merge phenotype, PGS, and CNV burden data
  df_temp <- left_join(pheno[, names(pheno) %in% c("IID", p)], cnv, by = "IID")
  colnames(df_temp) [2] <- "PHENO"

  # Fit linear regression - CNV burden
  fit_CNV <- lm(PHENO ~ CNV_GENES, data = df_temp)
  
  # Fill in the result table
  df[i, "EFFECT_CNV_GENES"] <- summary(fit_CNV)$coefficients[2,1]
  df[i, "SE_CNV_GENES"] <- summary(fit_CNV)$coefficients[2,2]
  df[i, "P_CNV_GENES"] <- summary(fit_CNV)$coefficients[2,4]
  
  df[i, "adjR2_CNV_GENES"] <- summary(fit_CNV)$adj.r.squared
  
  # Fit linear regression - DEL burden
  fit_DEL <- lm(PHENO ~ DEL_GENES, data = df_temp)
  
  # Fill in the result table
  df[i, "EFFECT_DEL_GENES"] <- summary(fit_DEL)$coefficients[2,1]
  df[i, "SE_DEL_GENES"] <- summary(fit_DEL)$coefficients[2,2]
  df[i, "P_DEL_GENES"] <- summary(fit_DEL)$coefficients[2,4]
  
  df[i, "adjR2_DEL_GENES"] <- summary(fit_DEL)$adj.r.squared
  
  # Fit linear regression - DUP burden
  fit_DUP <- lm(PHENO ~ DUP_GENES, data = df_temp)
  
  # Fill in the result table
  df[i, "EFFECT_DUP_GENES"] <- summary(fit_DUP)$coefficients[2,1]
  df[i, "SE_DUP_GENES"] <- summary(fit_DUP)$coefficients[2,2]
  df[i, "P_DUP_GENES"] <- summary(fit_DUP)$coefficients[2,4]
  
  df[i, "adjR2_DUP_GENES"] <- summary(fit_DUP)$adj.r.squared
  
}
rm(i, p, df_temp, fit_CNV, fit_DUP, fit_DEL)


########################################################
# STEP 3: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_phenotype/CNVBurden_model.txt", col.names = T, row.names = F, quote = F, sep = "\t")
