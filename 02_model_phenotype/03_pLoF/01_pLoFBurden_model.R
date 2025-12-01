# Model adjusted phenotypes according to pLoF burden for 43 assessed phenotypes

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)


########################################################
# STEP 1: Load data
########################################################

# Phenotypic variability explained by the PGS
# File with the following columns:  PHENO, cor_2 (squared correlation between phenotype and PGS)
cor2 <- as.data.frame(fread("CNV_PGS/data/phenotype_explained_by_PGS.txt"))

# Phenotype (covariate-corrected + INT; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing covariate-adjusted, inverse-normal transformed phenotype values 
pheno <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.txt"))
pheno <- pheno[pheno$IID %in% test_samples$IID, ]

# pLoF burden
# File with sample identifier (IID) as first column and the total pLoF burden as second column
lof <- as.data.frame(fread("CNV_PGS/data/LoF_burden.total.txt"))


########################################################
# STEP 2: Regression analysis
########################################################

# Create a dataframe to store regression results
df <- cor2[, c(1,5)]
rm(cor2)

# Loop over phenotypes
for (i in 1:nrow(df)) {
  
  # Define the signal
  p <- df[i, "PHENO"]

  # Merge phenotype and LoF burden data
  df_temp <- na.omit(left_join(pheno[, names(pheno) %in% c("IID", p)], lof, by = "IID"))
  colnames(df_temp)[2] <- c("PHENO")

  # Fit linear regression - pLoF burden
  fit <- lm(PHENO ~ LoF_total, data = df_temp)
  
  # Fill in the result table
  df[i, "EFFECT_LoF"] <- summary(fit)$coefficients[2,1]
  df[i, "SE_LoF"] <- summary(fit)$coefficients[2,2]
  df[i, "P_LoF"] <- summary(fit)$coefficients[2,4]
  
  df[i, "adjR2_LoF"] <- summary(fit)$adj.r.squared
  
}
rm(i, p, df_temp, fit)


########################################################
# STEP 3: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_phenotype/LoFBurden_model.txt", col.names = T, row.names = F, quote = F, sep = "\t")
