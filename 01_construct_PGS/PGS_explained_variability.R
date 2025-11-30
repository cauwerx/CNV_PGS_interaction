# Assess the proportion of phenotypic variability explained by the PGS

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(MBESS)


########################################################
# STEP 1: Load data
########################################################

# Testing samples; 
# File with a single column (IID) containing the sample identifier for all samples in the test set 
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# Phenotype (INT but not covariate-corrected; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing inverse-normal transformed phenotype values 
pheno_int <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_test_INT_All.tsv"))
pheno_int <- pheno_int[pheno_int$IID %in% test_samples$IID, ]

# Polygenic scores (PGSs; filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS values
pgs <- as.data.frame(fread("CNV_PGS/data/PGS_continuous.txt.gz"))
pgs <- pgs[pgs$IID %in% test_samples$IID, ]

# Covariates (filtered for testing samples IID) )
# File with sample identifier (IID) as first column, then one column per covariate, including age, age^2, sex, genotyping batch, and PC1-40
cov <- as.data.frame(fread("CNV_PGS/data/covariates.txt"))
cov <- cov[cov$IID %in% test_samples$IID, ]


########################################################
# STEP 2: Performance of PGS
########################################################

# Create a dataframe to store results
df <- data.frame(PHENO = names(pheno_int)[-1])

# Set parameters 
alpha <- 0.05

# Loop over phenotypes
for (p in df$PHENO) {
  
  print(paste0("Analyzing: ", p))

  # Subset data for a given phenotype across input data & merge
  df_pheno <- pheno_int[, c("IID", p)]; colnames(df_pheno)[2] <- "PHENO"
  df_pgs <- pgs[, c("IID", p)]; colnames(df_pgs)[2] <- "PGS"
  df_temp <- na.omit(left_join(df_pheno, df_pgs, by = "IID"))
  df_temp <- na.omit(left_join(df_temp, cov, by = "IID"))
  
  # Add number of samples, n
  n <- nrow(df_temp)
  df[which(df$PHENO == p), "N"] <- n
  
  ##########################################################
  # MODEL 1: Y_int ~ PGS + sex + age + age2 + batch + PC1-40
  fit1 <- lm(PHENO ~ . , data = df_temp[, -c(1)])
  
  # Fill in table with summary statistics
  df[which(df$PHENO == p), "fit1_BETA"] <- summary(fit1)$coefficients[2,1]
  df[which(df$PHENO == p), "fit1_SE"] <- summary(fit1)$coefficients[2,2]
  df[which(df$PHENO == p), "fit1_P"] <- summary(fit1)$coefficients[2,4]
  df[which(df$PHENO == p), "fit1_R2"] <- summary(fit1)$r.squared
  df[which(df$PHENO == p), "fit1_R2adj"] <- summary(fit1)$adj.r.squared
  
  # Add CI for R2
  k <- nrow(summary(fit1)$coefficients) - 1 # Remove the intercept, which is not a predictor
  df[which(df$PHENO == p), "fit1_R2_lower"] <- ci.R2(R2 = df[which(df$PHENO == p), "fit1_R2"], N = n, K = k, conf.level = 0.95)[[1]]
  df[which(df$PHENO == p), "fit1_R2_upper"] <- ci.R2(R2 = df[which(df$PHENO == p), "fit1_R2"], N = n, K = k, conf.level = 0.95)[[3]]

  
  ##########################################################
  # MODEL 2: Y_int ~ sex + age + age2 + batch + PC1-40 (NO PGS)
  fit2 <- lm(PHENO ~ . , data = df_temp[, -c(1,3)])
  
  # Fill in table with summary statistics
  df[which(df$PHENO == p), "fit2_R2"] <- summary(fit2)$r.squared
  df[which(df$PHENO == p), "fit2_R2adj"] <- summary(fit2)$adj.r.squared
  
  # Add CI for R2
  k <- ncol(df_temp[, -c(1,3)]) - 1 # Remove the intercept, which is not a predictor
  df[which(df$PHENO == p), "fit2_R2_lower"] <- ci.R2(R2 = df[which(df$PHENO == p), "fit2_R2"], N = n, K = k, conf.level = 0.95)[[1]]
  df[which(df$PHENO == p), "fit2_R2_upper"] <- ci.R2(R2 = df[which(df$PHENO == p), "fit2_R2"], N = n, K = k, conf.level = 0.95)[[3]]
  
  
  ##########################################################
  # DELTA: Difference in R^2 between model 1 and 2
  
  # Calculate Delta R^2, the incremental improvement in explained phenotypic variability 
  df[which(df$PHENO == p), "DELTA_R2_fit1_fit2"] <- df[which(df$PHENO == p), "fit1_R2"] - df[which(df$PHENO == p), "fit2_R2"]
  
  # Add CI for DELTA
  k <- 1 
  df[which(df$PHENO == p), "DELTA_R2_lower"] <- ci.R2(R2 = df[which(df$PHENO == p), "DELTA_R2_fit1_fit2"], N = n, K = k, conf.level = 0.95)[[1]]
  df[which(df$PHENO == p), "DELTA_R2_upper"] <- ci.R2(R2 = df[which(df$PHENO == p), "DELTA_R2_fit1_fit2"], N = n, K = k, conf.level = 0.95)[[3]]

}
rm(p, df_pheno, df_pgs, df_temp, fit1, fit2, n, k)


# Save file
fwrite(df, "phenotype_explained_by_PGS.txt", col.names = T, row.names = F, quote = F, sep = "\t")


########################################################
# STEP 3: Assess Delta R^2
########################################################

# Calculate median Delta
median_delta <- median(df$DELTA_R2_fit1_fit2)

# Calculate IQR Delta
iqr_delta <- IQR(df$DELTA_R2_fit1_fit2)

# Identify Delta outliers
df[which(df$DELTA_R2_fit1_fit2 > 0.2), )]
df[which(df$DELTA_R2_fit1_fit2 < 0.05), )]