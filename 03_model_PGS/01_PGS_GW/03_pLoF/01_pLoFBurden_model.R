# Model genome-wide PGS_gw according to the pLoF burden for traits which are significantly impacted by the pLoF burden

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

# Phenotypic variability explained by the PGS
# File with the following columns:  PHENO, cor_2 (squared correlation between phenotype and PGS)
cor2 <- as.data.frame(fread("CNV_PGS/data/phenotype_explained_by_PGS.txt"))

# Polygenic scores (PGSs; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous.txt.gz"))
pgs <- pgs[pgs$IID %in% test_samples$IID, ]

# pLoF burden
# File with sample identifier (IID) as first column and the total pLoF burden as second column
lof <- as.data.frame(fread("CNV_PGS/data/LoF_burden.total.txt"))
lof <- lof[lof$IID %in% test_samples$IID, ]

# Load LoF burden effect on phenotype; Output from model_phenotype/pLoFBurden_model.R 
df_lof_on_pheno <- as.data.frame(fread("CNV_PGS/data/model_phenotype/LoFBurden_model.txt"))


########################################################
# STEP 2: Regression analysis
########################################################

# Select traits to test
pheno_lof_burden <- df_lof_on_pheno[which(df_lof_on_pheno$P_LoF <= 0.05/43), "PHENO"] # 15

# Create a dataframe to store results
df <- df_lof_on_pheno[, c(1:5)]

# Loop over selected phenotypes - impacted by the pLoF burden
for (p in pheno_lof_burden) {
  
  # Merge PGS and pLoF burden data
  df_temp <- na.omit(left_join(pgs[, names(pgs) %in% c("IID", p)], lof, by = "IID"))
  colnames(df_temp) [2] <- "PGS"
  
  # Fit linear regression
  fit <- lm(PGS ~ LoF_total, data = df_temp[, -c(1)])
  
  # Fill table with results
  df[which(df$PHENO == p), "EFFECT_LoF_on_PGS_GW"] <- summary(fit)$coefficients[2,1]
  df[which(df$PHENO == p), "SE_LoF_on_PGS_GW"] <- summary(fit)$coefficients[2,2]
  df[which(df$PHENO == p), "P_LoF_on_PGS_GW"] <- summary(fit)$coefficients[2,4]

}
rm(p, df_temp, fit)


########################################################
# STEP 3: Compare to CNV burden effect on phenotype
########################################################

# Sign concordance (i.e., if CNV burde increases a trait, the CNV burden will tend to increase the PGS for that trait & vice versa)

# Among nominally significant 
df_nom <- df[which(df$P_LoF_on_PGS_GW <= 0.05), ]
nrow(df_nom[which((df_nom$EFFECT_LoF > 0 & df_nom$EFFECT_LoF_on_PGS_GW > 0) | (df_nom$EFFECT_LoF < 0 & df_nom$EFFECT_LoF_on_PGS_GW < 0)), ])
# All 1 are directionally concordant
binom.test(1,1,0.5, alternative = "greater")


########################################################
# STEP 4: Save data
########################################################

fwrite(df, "data/model_PGS/LoFBurden_model_PGS.txt", col.names = T, row.names = F, quote = F, sep = "\t")
