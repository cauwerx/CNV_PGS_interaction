# Model raw phenotypes according to CNV status, PGS, and their interaction for 119 CNV-trait pairs

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(tidyr)


########################################################
# STEP 1: Load data
########################################################

# Testing samples
# File with a single column (IID) containing the sample identifier for all samples in the test set
test_samples <- as.data.frame(fread("CNV_PGS/data/test_IDs.txt"))

# CNV signals
# File with the following columns: PHENO, CHR, CNVR_START, CNVR_STOP, TOP_MODEL, CB (cytogenic band)
cnv_signals <- as.data.frame(fread("CNV_PGS/data/cnv_signals.txt"))

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
rm(cnv_signals)

# Loop over phenotypes
for (i in 1:nrow(df)) {
  
  # Define the signal
  p <- df[i, "PHENO"]
  chr <- df[i, "CHR"]
  pos <-  df[i, "TOP_POS"]
  model <- df[i, "TOP_MODEL"]; if(model == "M-DUP") {model <- "M"}
  print(paste0("Analyzing ", p," for chr", chr, ":", pos, " CNVs (", model, ")"))
  
  # Identify samples carrying relevant CNVs (overlapping the lead probe)
  df_cnvs <- cnvs[which(cnvs$Chromosome == chr & cnvs$Start_Position_bp <= pos & cnvs$End_Position_bp >= pos), ]
  
  # Identify high-confidence deletion and duplication carriers
  del_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  dup_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
  
  # Merge raw phenotype & PGS data
  df_temp <- left_join(pheno[, names(pheno) %in% c("IID", p)], pgs[, names(pgs) %in% c("IID", p)], by = "IID")
  colnames(df_temp) <- c("IID", "PHENO", "PGS")
  
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
  df_temp <- left_join(df_temp, cov, by = "IID")
  
  # Fit linear regression
  fit <- lm(PHENO ~ . + CNV:PGS, data = df_temp[, -c(1)])
  
  # Fill in the result table with sample number and regression coefficients
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
  
  df[i, "EFFECT_CNV"] <- summary(fit)$coefficients[3,1]
  df[i, "SE_CNV"] <- summary(fit)$coefficients[3,2]
  df[i, "P_CNV"] <- summary(fit)$coefficients[3,4]
  
  df[i, "EFFECT_PGS"] <- summary(fit)$coefficients[2,1]
  df[i, "SE_PGS"] <- summary(fit)$coefficients[2,2]
  df[i, "P_PGS"] <- summary(fit)$coefficients[2,4]
  
  df[i, "EFFECT_CNVxPGS"] <- summary(fit)$coefficients[152,1]
  df[i, "SE_CNVxPGS"] <- summary(fit)$coefficients[152,2]
  df[i, "P_CNVxPGS"] <- summary(fit)$coefficients[152,4]
  
  df[i, "adjR2"] <- summary(fit)$adj.r.squared
  
}
rm(i, p, chr, pos, model, del_carriers, dup_carriers, df_temp)


########################################################
# STEP 3: Save data
########################################################

fwrite(df, "CNV_PGS/data/model_phenotype/CNV_PGS_interaction_model_RawScale.txt", col.names = T, row.names = F, quote = F, sep = "\t")


########################################################################################################################################################################

########################################################
# STEP 4: Follow-up: significant interactions
########################################################
# Test whether the interaction survives when using PGS_trans_250kb instead of PGS_gw

# Polygenic scores (PGSs trans 250kb (i.e., excluding SNVs ± 250kb around the CNV region); filtered for testing samples IID) 
# File with sample identifier (IID) as first column, then one column per phenotype, containing PGS trans 250 kb values
pgs_trans <- as.data.frame(fread("CNV_PGS/data/PGS_continuous_out_250kb.txt.gz"))
pgs_trans <- pgs_trans[pgs_trans$IID %in% test_samples$IID, ]


#########################
# Analyze 22q11.23 & GGT

# Define the signal
p <- "GGT"
chr <- 22
pos <-  24836024
start <- 23688345
end <- 24990213
model <- "M"

# Identify samples carrying relevant CNVs (overlapping the lead probe)
df_cnvs <- cnvs[which(cnvs$Chromosome == chr & cnvs$Start_Position_bp <= pos & cnvs$End_Position_bp >= pos), ]

# Identify high-confidence deletion and duplication carriers
del_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
dup_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]

# Merge phenotype, PGS trans 250kb, CNV carrier status, and covariate information
df_temp <- left_join(pheno[, names(pheno) %in% c("IID", p)], pgs_trans[, names(pgs_trans) %in% c("IID", paste0(p, "_chr", chr, "_", start-250000, "_", end+250000))], by = "IID")
colnames(df_temp) <- c("IID", "PHENO", "PGS_trans")
df_temp$CNV <- 0
df_temp[df_temp$IID %in% del_carriers, "CNV"] <- -1
df_temp[df_temp$IID %in% dup_carriers, "CNV"] <- 1
df_temp <- na.omit(left_join(df_temp, cov, by = "IID"))

# Fit linear regression & print results
fit <- lm(PHENO ~ . + CNV:PGS_trans, data = df_temp[, -c(1)])
summary(fit)
rm(p, chr, pos, start, end, model, df_cnvs, del_carriers, dup_carriers, df_temp, fit)

########################
# Analyze 22q11.23 & GS

# Define the signal
p <- "GS"
chr <- 22
pos <-  24787635
start <- 23688345
end <- 24985853
model <- "DUP"

# Identify samples carrying relevant CNVs (overlapping the lead probe)
df_cnvs <- cnvs[which(cnvs$Chromosome == chr & cnvs$Start_Position_bp <= pos & cnvs$End_Position_bp >= pos), ]

# Identify high-confidence deletion and duplication carriers
del_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 1 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]
dup_carriers <- df_cnvs[which(df_cnvs$Copy_Number == 3 & abs(df_cnvs$Quality_Score) >= 0.5), "IID"]

# Merge phenotype, PGS trans 250kb, CNV carrier status, and covariate information
df_temp <- left_join(pheno[, names(pheno) %in% c("IID", p)], pgs_trans[, names(pgs_trans) %in% c("IID", paste0(p, "_chr", chr, "_", start-250000, "_", end+250000))], by = "IID")
colnames(df_temp) <- c("IID", "PHENO", "PGS_trans")
df_temp$CNV <- 0
df_temp[df_temp$IID %in% dup_carriers, "CNV"] <- 1
df_temp <- na.omit(left_join(df_temp, cov, by = "IID"))

# Fit linear regression & print results
fit <- lm(PHENO ~ . + CNV:PGS_trans, data = df_temp[, -c(1)])
summary(fit)
rm(p, chr, pos, start, end, model, df_cnvs, del_carriers, dup_carriers, df_temp, fit)
