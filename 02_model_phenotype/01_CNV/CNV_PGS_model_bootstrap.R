# Model adjusted phenotypes according to CNV status VS according to CNV status and PGS for 119 CNV-trait pairs
# Correlation between the CNV effects from the two fitted models is estimated through bootstrap


########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)
library(modelr)
library(purrr)


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

# Create a dataframe to store resgression results
df <- right_join(cor2, cnv_signals, by = "PHENO")
rm(cor2)

# Set parameters
set.seed(678)
n_batch <- 5
batch_size <- 1000

# Loop over CNV-trait pairs
for (i in 1:nrow(df)) {
  
  ####################################
  # Set up the data
  ####################################
  
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
  
  # Merge pheno and PGS data
  df_temp <- left_join(pheno[, names(pheno) %in% c("IID", p)], pgs[, names(pgs) %in% c("IID", p)], by = "IID")
  colnames(df_temp) <- c("IID", "PHENO", "PGS")

  #  Add CNV carrier status for the relevant model
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
  

  ####################################
  # Fit models & store data
  ####################################
  
  ############################################
  # Model 1: CNV only
  
  # Fit linear regression 
  fit1 <- lm(PHENO ~ CNV, data = df_temp)
  
  # Fill in the result table
  df[i, "M1_EFFECT_CNV"] <- summary(fit1)$coefficients[2,1]
  df[i, "M1_SE_CNV"] <- summary(fit1)$coefficients[2,2]
  df[i, "M1_P_CNV"] <- summary(fit1)$coefficients[2,4]

  
  ############################################
  # M2: CNV and PGS model
  
  # Fit linear regression 
  fit2 <- lm(PHENO ~ CNV + PGS, data = df_temp)
  
  # Fill in the result table
  df[i, "M2_EFFECT_CNV"] <- summary(fit2)$coefficients[2,1]
  df[i, "M2_SE_CNV"] <- summary(fit2)$coefficients[2,2]
  df[i, "M2_P_CNV"] <- summary(fit2)$coefficients[2,4]
  
  df[i, "M2_EFFECT_PGS"] <- summary(fit2)$coefficients[3,1]
  df[i, "M2_SE_PGS"] <- summary(fit2)$coefficients[3,2]
  df[i, "M2_P_PGS"] <- summary(fit2)$coefficients[3,4]
 

  ####################################
  # Bootstrap
  # https://svmiller.com/blog/2020/03/bootstrap-standard-errors-in-r/
  ####################################
  
  # Bootstrap resample the dataset, run lm(), summarize, extract effect size estimates, store in a single table -> M1 (CNV-only)
  df_boot1_tidy <- data.frame()
  for (b in 1:n_batch) { # Split in batches to avoid running out of RAM
    # Bootstrap round
    df_boot1 <- df_temp %>%
                bootstrap(batch_size) %>%
                mutate(lm = map(strap, ~lm(PHENO ~ CNV, data = .)),
                       tidy = map(lm, broom::tidy)) %>%
                pull(tidy) %>%
                map2_df(seq(1, batch_size), ~mutate(.x, resample = .y))
    # Save beta estimates
    df_boot1 <- as.data.frame(df_boot1[which(df_boot1$term == "CNV"), c(1:2)])
    # Merge to summary table
    df_boot1_tidy <- rbind(df_boot1_tidy, df_boot1)
    # Clean
    rm(df_boot1)
  }
  
  # Bootstrap resample the dataset, run lm(), summarize, extract effect size estimates, store in a single table -> M2 (CNV & PGS)
  df_boot2_tidy <- data.frame()
  for (b in 1:n_batch) { # Split in batches to avoid running out of RAM
    # Bootstrap round
    df_boot2 <- df_temp %>%
      bootstrap(batch_size) %>%
      mutate(lm = map(strap, ~lm(PHENO ~ CNV + PGS, data = .)),
             tidy = map(lm, broom::tidy)) %>%
      pull(tidy) %>%
      map2_df(seq(1, batch_size), ~mutate(.x, resample = .y))
    # Save beta estimates
    df_boot2 <- as.data.frame(df_boot2[which(df_boot2$term == "CNV"), c(1:2)])
    # Merge to summary table
    df_boot2_tidy <- rbind(df_boot2_tidy, df_boot2)
    # Clean
    rm(df_boot2)
  }

  # Get the correlation between boostrap estimates & save to results
  df[i, "bootstrap_r"] <- cor(df_boot1_tidy$estimate, df_boot2_tidy$estimate, method = "pearson")
  
}
rm(i, p, chr, pos, model, del_carriers, dup_carriers, df_temp, fit1, fit2, df_boot, df_boot1, df_boot2, df_boot1_tidy, df_boot2_tidy)


########################################################
# STEP 3: Compare effect sizes
########################################################

# Add significant differences in CNV effects between M1 (CNV-only predictor) and M2 (CNV & PGS predictors)
df$T_diff_effect_CNV_M1_M2 <- (df$M1_EFFECT_CNV - df$M2_EFFECT_CNV)/sqrt(df$M1_SE_CNV^2 + df$M2_SE_CNV^2 - (2*df$M1_SE_CNV*df$M2_SE_CNV*df$bootstrap_r))
df$P_diff_effect_CNV_M1_M2 <- 2*pnorm(-abs(df$T_diff_effect_CNV_M1_M2), mean = 0, sd = 1)


########################################################
# STEP 4: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_phenotype/CNV_vs_CNV_PGS_model_bootstrap.txt", col.names = T, row.names = F, quote = F, sep = "\t")