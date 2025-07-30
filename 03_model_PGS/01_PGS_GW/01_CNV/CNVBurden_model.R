# Model genome-wide PGS_gw according to CNV burden for traits which are significantly impacted by the CNV burden

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

# Polygenic scores (PGSs; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous.txt.gz"))
pgs <- pgs[pgs$IID %in% test_samples$IID, ]

# CNV burden data (filtered for testing samples IID)
# CNV burdens for the UKBB (as described in Auwerx et al., 2022)
cnv_burden <- as.data.frame(fread("CNV_PGS/data/CNV_burden.txt.gz", select = c(1,10), col.names = c("IID", "CNV_GENES")))
dup_burden <- as.data.frame(fread("CNV_PGS/data/DUP_burden.txt.gz", select = c(1,10), col.names = c("IID", "DUP_GENES")))
del_burden <- as.data.frame(fread("CNV_PGS/data/DEL_burden.txt.gz", select = c(1,10), col.names = c("IID", "DEL_GENES")))
cnv <- left_join(cnv_burden, dup_burden, by = "IID")
cnv <- left_join(cnv, del_burden, by = "IID")
cnv <- cnv[cnv$IID %in% test_samples$IID, ]
rm(cnv_burden, dup_burden, del_burden)

# Load CNV burden effect on phenotype; Output from model_phenotype/CNVBurden_model.R
df_cnv_on_pheno <- as.data.frame(fread("CNV_PGS/data/model_phenotype/CNVBurden_model.txt"))


########################################################
# STEP 2: Regression analysis
########################################################

# Select traits to test
pheno_cnv_burden <- df_cnv_on_pheno[which(df_cnv_on_pheno$P_CNV_GENES <= 0.05/43), "PHENO"]

# Create a dataframe to store results
df <- df_cnv_on_pheno[, c(1:5)]
colnames(df) <- c("PHENO", "cor_2", "EFFECT_CNV_BURDEN", "SE_CNV_BURDEN", "P_CNV_BURDEN")
df$CNV_TYPE <- "CNV"

# Loop over selected phenotypes - impacted by the CNV burden
for (p in pheno_cnv_burden) {
  
  # Define the signal
  type <- "CNV"

  # Merge PGS and CNV burden data
  df_temp <- na.omit(left_join(pgs[, names(pgs) %in% c("IID", p)], cnv, by = "IID"))
  colnames(df_temp) [2] <- "PGS"
  
  # Fit linear regression
  fit <- lm(PGS ~ CNV_GENES, data = df_temp)
  
  # Fill table with results
  df[which(df$PHENO == p & df$CNV_TYPE == type), "EFFECT_CNV_BURDEN_on_PGS_GW"] <- summary(fit)$coefficients[2,1]
  df[which(df$PHENO == p & df$CNV_TYPE == type), "SE_CNV_BURDEN_on_PGS_GW"] <- summary(fit)$coefficients[2,2]
  df[which(df$PHENO == p & df$CNV_TYPE == type), "P_CNV_BURDEN_on_PGS_GW"] <- summary(fit)$coefficients[2,4]
  
}
rm(type, p, df_temp, fit)


########################################################
# STEP 3: Compare to CNV burden effect on phenotype
########################################################

# Sign concordance (i.e., if CNV burde increases a trait, the CNV burden will tend to increase the PGS for that trait & vice versa)

# Among all CNV-trait pairs
nrow(df[which((df$EFFECT_CNV_BURDEN > 0 & df$EFFECT_CNV_BURDEN_on_PGS_GW > 0) | (df$EFFECT_CNV_BURDEN < 0 & df$EFFECT_CNV_BURDEN_on_PGS_GW < 0)), ])
# 18 are directionally concordant
binom.test(18,21,0.5, alternative = "greater")

# Among nominally significant 
df_nom <- df[which(df$P_CNV_BURDEN_on_PGS_GW <= 0.05), ]
nrow(df_nom[which((df_nom$EFFECT_CNV_BURDEN > 0 & df_nom$EFFECT_CNV_BURDEN_on_PGS_GW > 0) | (df_nom$EFFECT_CNV_BURDEN < 0 & df_nom$EFFECT_CNV_BURDEN_on_PGS_GW < 0)), ])
# All 4 are directionally concordant
binom.test(4,4,0.5, alternative = "greater")


########################################################
# STEP 4: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_PGS/CNVBurden_model_PGS.txt", col.names = T, row.names = F, quote = F, sep = "\t")