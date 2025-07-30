# Model raw phenotypes according to PGS among individuals carrying a relevant CNV for 119 CNV-trait pairs


########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)
library(readxl)
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

# Phenotype (raw; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then one column per phenotype, containing raw phenotype values 
pheno <- as.data.frame(fread("CNV_PGS/data/pheno_continuous_WB_raw.txt"))
pheno <- pheno[pheno$IID %in% test_samples$IID, ]

# Covariates (raw; filtered for testing samples IID)
# File with sample identifier (IID) as first column, then age, age2^2, sex, array (i.e., genotyping array), batch (i.e., genotyping batch), principal components 1-40
cov <- as.data.frame(fread("CNV_PGS/data/covariates.txt"))
cov <- cov[cov$IID %in% test_samples$IID, ]

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
df <- cnv_signals
colnames(df)[9:11] <- c("EFFECT_CNV", "SE_CNV", "P_CNV")

# Add phenotypic variance explained by PGS
df <- right_join(cor2, df, by = "PHENO")
rm(cor2)

# Set parameters
set.seed(678)
n_bootstrap <- 5000

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
  scaling_factor <- df[i, "cor_2"]/sqrt(df[i, "var_PGS"])
  print(paste0("Analyzing ", p," for chr", chr, ":", pos, " CNVs (", model, ")"))
  
  # Identify samples carrying relevant CNVs (overlaping lead signal)
  df_cnvs <- cnvs[which(cnvs$Chromosome == chr & cnvs$Start_Position_bp <= pos & cnvs$End_Position_bp >= pos), ]
  
  # Identify high confidence deletion and duplication carriers
  del_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  dup_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  
  # Create a temporary dataframe
  df_temp <- data.frame(IID = test_samples$IID)

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
    df_temp[df_temp$IID %in% del_carriers, "CNV"] <- 1
    df_temp[df_temp$IID %in% dup_carriers, "CNV"] <- 1
  }

  # Add phenotype values and PGS
  df_temp <- left_join(df_temp, pheno[, names(pheno) %in% c("IID", p)], by = "IID")
  colnames(df_temp)[3] <- "PHENO"
  df_temp <- left_join(df_temp, pgs[, names(pgs) %in% c("IID", p)], by = "IID")
  colnames(df_temp)[4] <- "PGS"

  # Add covariates (age, age^2, sex, PC1-10) & factorize (sex)
  df_temp <- left_join(df_temp, cov[, -c(5:6, 17:46)], by = "IID")
  df_temp$sex <- factor(df_temp$sex)
  
  # Get the population SD
  pop_sd <- sd(df_temp$PHENO, na.rm = T)
  
  # Subset CNV carriers from control copy-neutral individuals
  df_carrier <- na.omit(df_temp[which(df_temp$CNV == 1), !names(df_temp) %in% c("IID", "CNV")])

  
  #############################
  # BOOTSTRAP: carriers
  # https://svmiller.com/blog/2020/03/bootstrap-standard-errors-in-r/
  #############################
  
  # Fit a linear regression 
  fit_carrier <- lm(PHENO ~ ., data = df_carrier)
  
  # Bootstrap resample the dataset
  boot_carrier <- df_carrier %>%
                  bootstrap(n_bootstrap)
  
  # Run lm() on each bootstrap resample and store results
  boot_carrier <- boot_carrier %>% 
                  mutate(lm = map(strap, ~lm(PHENO ~ ., data = .)), tidy = map(lm, broom::tidy))

  # Summarize results in a single table
  boot_carrier_tidy <- boot_carrier %>%
                        pull(tidy) %>%
                        map2_df(., seq(1, n_bootstrap), ~mutate(.x, resample = .y)) 
  
  # Calculate the bootstrap SE, which corresponds to the SD of the coefficient estimate in the model
  BSE_carrier <- boot_carrier_tidy %>%
                 group_by(term) %>%
                 summarize(bse = sd(estimate))
  
  
  # Get the p-value for the coeffcients
  fit_carrier_summary <- broom::tidy(fit_carrier) %>%
                          mutate(category = "OLS Standard Errors") %>%
                          bind_rows(., broom::tidy(fit_carrier) %>% select(-std.error) %>% left_join(., BSE_carrier %>% rename(std.error = bse))) %>%
                          mutate(category = ifelse(is.na(category), "Bootstrapped Standard Errors", category),
                                 tstat = estimate/std.error,
                                 pval = 1.96*pt(-abs(tstat), df = fit_carrier$df), # dfs from fit_carrier
                                 L95 = estimate - 1.96*std.error,
                                 U95 = estimate + 1.96*std.error)
  
  # Fill in the result table
  df[i, "N_carrier"] <- nrow(df_carrier)
  df[i, "population_SD"] <- pop_sd
  df[i, "EFFECT_PGS_carrier"] <- fit_carrier_summary[which(fit_carrier_summary$category == "Bootstrapped Standard Errors" & fit_carrier_summary$term == "PGS"), "estimate"]
  df[i, "SE_PGS_carrier"] <- fit_carrier_summary[which(fit_carrier_summary$category == "Bootstrapped Standard Errors" & fit_carrier_summary$term == "PGS"), "std.error"]
  df[i, "P_PGS_carrier"] <- fit_carrier_summary[which(fit_carrier_summary$category == "Bootstrapped Standard Errors" & fit_carrier_summary$term == "PGS"), "p.value"]
  df[i, "BootSE_PGS_carrier"] <- BSE_carrier[which(BSE_carrier$term == "PGS"), "bse"]
  df[i, "BootP_PGS_carrier"] <- fit_carrier_summary[which(fit_carrier_summary$category == "Bootstrapped Standard Errors" & fit_carrier_summary$term == "PGS"), "pval"]
  df[i, "BootL95_PGS_carrier"] <- fit_carrier_summary[which(fit_carrier_summary$category == "Bootstrapped Standard Errors" & fit_carrier_summary$term == "PGS"), "L95"]
  df[i, "BootU95_PGS_carrier"] <- fit_carrier_summary[which(fit_carrier_summary$category == "Bootstrapped Standard Errors" & fit_carrier_summary$term == "PGS"), "U95"]
  
  # For Figure S2, plotting effects on the SD scale
  df[i, "EFFECT_PGS_SD"] <- df[i, "EFFECT_PGS_carrier"]/df[i, "population_SD"]
  df[i, "L95_PGS_SD"] <- quantile(boot_carrier_tidy$estimate[which(boot_carrier_tidy$term == "PGS")], 0.025)/df[i, "population_SD"]
  df[i, "U95_PGS_SD"] <- quantile(boot_carrier_tidy$estimate[which(boot_carrier_tidy$term == "PGS")], 0.975)/df[i, "population_SD"]
  
  rm(df_carrier, fit_carrier, boot_carrier, boot_carrier_tidy, BSE_carrier, fit_carrier_summary)
                 
}
rm(i, p, chr, pos, model, scaling_factor, df_cnvs, df_temp, fit)


########################################################
# STEP 4: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_phenotype/bootstrapSE_effect_PGS_carrier_vs_non_carrier_best_CNV_model.txt", col.names = T, row.names = F, quote = F, sep = "\t")


########################################################
# STEP 4: Assess effect of sample size and PGS performance on likelihood to be significant
########################################################

# Identify Bonferroni-significant signals
df$SIG <- 0
df[which(df$BootP_PGS_carrier <= 0.05/119), "SIG"] <- 1
df$SIG <- factor(df$SIG, levels = c(0,1))

# Fit a linear regression model
fit_carriers_perf <- glm(SIG ~ N_carrier + cor_2, data = df, family = "binomial")
summary(fit_carriers_perf)
