# Assess the proportion of phenotypic variability explained by the PGS

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)


########################################################
# STEP 1: Load data
########################################################

# Testing samples; 
# File with a single column (IID) containing the sample identifier for all samples in the test set. 
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# Phenotype (covariate-corrected + INT; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing covariate-adjusted, inverse-normal transformed phenotype values 
pheno <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.txt"))
pheno <- pheno[pheno$IID %in% test_samples$IID, ]

# Polygenic scores (PGSs; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous.txt.gz"))
pgs <- pgs[pgs$IID %in% test_samples$IID, ]


########################################################
# STEP 2: Performance of PGS
########################################################

# Create a dataframe to store results
df <- data.frame(PHENO = names(pheno)[-1])

# Loop over phenotypes
for (p in df$PHENO) {
  
  # Subset the data
  df_pheno <- pheno[, c("IID", p)]; colnames(df_pheno)[2] <- "PHENO"
  df_pgs <- pgs[, c("IID", p)]; colnames(df_pgs)[2] <- "PGS"
  df_temp <- na.omit(left_join(df_pheno, df_pgs, by = "IID")) 
  df[which(df$PHENO == p), "cor_2"] <- cor(df_temp[,2], df_temp[,3])^2
  df[which(df$PHENO == p), "var_PGS"] <- var(df_temp[,"PGS"])
  print(paste0("Correlation between ", p, " and its PGS in ", nrow(df_temp), " individuals: ", round(cor(df_temp[,2], df_temp[,3]),3)))

}
rm(p, df_pheno, df_pgs, df_temp)

# Order by explained variance
df <- df[order(df$cor_2), ]


########################################################
# STEP 3: Save data
########################################################

# Save file
fwrite(df, "CNV_PGS/data/phenotype_explained_by_PGS.txt", col.names = T, row.names = F, quote = F, sep = "\t")