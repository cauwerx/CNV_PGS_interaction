############################################################################################################
# Libraries
############################################################################################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)

############################################################################################################
# STEP 1: Load Main Data
############################################################################################################

DIR="/home_dir/"
OUT_DIR="/project/data/phenotypes/"

# NOTE! use test IDs only
pheno_all <- as.data.frame(fread(paste0(DIR, "/project/data/phenotypes/pheno_continuous_test_raw.tsv"), header = T))

pheno_M <- subset(pheno_all, sex==1)
pheno_M <- pheno_M[, !names(pheno_M) %in% c("menarche", "menopause","birth_weight_first_child")]
pheno_M$WHRadjBMI <- residuals(lm(pheno_M$WHR ~ pheno_M$BMI, na.action = na.exclude))
print(paste0("Dimensions of phenotype table for white British males: ", ncol(pheno_M)-2, " pheno x ", nrow(pheno_M), " males"))

pheno_F <- subset(pheno_all, sex==2)
pheno_F <- pheno_F[, !names(pheno_F) %in% c("facial_hair", "balding")]
pheno_F$WHRadjBMI <- residuals(lm(pheno_F$WHR ~ pheno_F$BMI, na.action = na.exclude))
print(paste0("Dimensions of phenotype table for white British females: ", ncol(pheno_F)-2, " pheno x ", nrow(pheno_F), " females"))

# remove sex from pheno files (as already in cov)
pheno_all <- subset(pheno_all, select = -sex)
pheno_M <- subset(pheno_M, select = -sex)
pheno_F <- subset(pheno_F, select = -sex)

############################################################################################################
# STEP 2: Inverse Normal Transformation 
############################################################################################################

#  All
pheno_int_all <- data.frame(IID = pheno_all$IID)
for (p in 2:ncol(pheno_all)) {
        pheno <- colnames(pheno_all)[p]
        print(paste0("INT - All: ", pheno))
        pheno_int_all[, pheno] <- qnorm((rank(pheno_all[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_all[, p])))
}; rm (p, pheno)

fwrite(pheno_int_all, paste0(OUT_DIR, "/pheno_continuous_test_INT_All.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Males
pheno_int_M <- data.frame(IID = pheno_M$IID)
for (p in 2:ncol(pheno_M)) {
        pheno <- colnames(pheno_M)[p]
        print(paste0("INT - Males: ", pheno))
        pheno_int_M[, pheno] <- qnorm((rank(pheno_M[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_M[, p])))
}; rm (p, pheno)

fwrite(pheno_int_M, paste0(OUT_DIR, "/pheno_continuous_test_INT_M.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Females
pheno_int_F <- data.frame(IID = pheno_F$IID)
for (p in 2:ncol(pheno_F)) {
        pheno <- colnames(pheno_F)[p]
        print(paste0("INT - Females: ", pheno))
        pheno_int_F[, pheno] <- qnorm((rank(pheno_F[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_F[, p])))
}; rm (p, pheno)

fwrite(pheno_int_F, paste0(OUT_DIR, "/pheno_continuous_test_INT_F.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

###########################################################################################################
# STEP 3: Correct for covariates
############################################################################################################

## Read in covariates
cov <- as.data.frame(fread(paste0(DIR, "/data/covariates.txt"), header = T, select = c(1:4, 6:46)))
colnames(cov)[1] <- "IID"

# Add home location (1 km resolution)
location <- as.data.frame(fread(paste0(DIR, "/project/data/phenotypes/location_coord.tsv"), header = T))
colnames(location) <- c("IID", "east_1km", "north_1km")

cov <- merge(cov, location, by='IID', all.x=TRUE)

# Select subsets
# all individuals
cov_all <- cov[cov$IID %in% pheno_all$IID, ] 
print(paste0("Covariates all: ", ncol(cov_all)-1))
print(colnames(cov_all[-1]))
print(paste0("# individuals:", nrow(cov_all)))
# males
cov_M <- cov[cov$IID %in% pheno_M$IID, -c(4)] # remove sex from cov
print(paste0("Covariates males: ", ncol(cov_M)-1))
print(colnames(cov_M[-1]))
print(paste0("# males:", nrow(cov_M)))
# females
cov_F <- cov[cov$IID %in% pheno_F$IID, -c(4)] # remove sex from cov
print(paste0("Covariates females: ", ncol(cov_F)-1))
print(colnames(cov_F[-1]))
print(paste0("# females:", nrow(cov_F)))

# PCs
pcs <- paste0("PC", 1:40)
pcs_sum <- paste(pcs, collapse = " + ")

# Regress out covariates - All
pheno_int_cor_all <- data.frame(IID = pheno_int_all$IID)
for (p in 2:ncol(pheno_int_all)) {
        # Define phenotype
        pheno <- colnames(pheno_int_all)[p]
        print(paste0("Covariates - All: ", pheno))
        # Create a temporary dataframe: join all cov with the p phenotype
        temp <- left_join(pheno_int_all[, c(1,p)], cov_all, by = "IID") # NOTE! order of IIDs is always of the left df (also for right_join)
	    # change name of p phenotype (col 2 in temp) to "pheno"
        colnames(temp)[2] <- "pheno"
        # Regress out covariates, remove ID col
	    formula <- as.formula(paste0("pheno ~ age + age2 + sex + batch + east_1km + north_1km + ", pcs_sum))
        pheno_int_cor_all[, pheno] <- residuals(lm(formula , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_all, paste0(OUT_DIR, "/pheno_continuous_test_INT_age_age2_sex_batch_location_PCs_All.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Regress out covariates - Males
pheno_int_cor_M <- data.frame(IID = pheno_int_M$IID)
for (p in 2:ncol(pheno_int_M)) {
        # Define phenotype
        pheno <- colnames(pheno_int_M)[p]
        print(paste0("Covariates - Males: ", pheno))
        # Create a temporary dataframe: join all cov with the p phenotype
        temp <- left_join(pheno_int_M[, c(1,p)], cov_M, by = "IID")
	    # change name of p phenotype (col 2 in temp) to "pheno"
        colnames(temp)[2] <- "pheno"
        # Regress out covariates
	    formula <- as.formula(paste0("pheno ~ age + age2 + batch + east_1km + north_1km + ", pcs_sum))
        pheno_int_cor_M[, pheno] <- residuals(lm(formula , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_M, paste0(OUT_DIR, "/pheno_continuous_test_INT_age_age2_batch_location_PCs_M.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Regress out covariates - Females
pheno_int_cor_F <- data.frame(IID = pheno_int_F$IID)
for (p in 2:ncol(pheno_int_F)) {
        # Define phenotype
        pheno <- colnames(pheno_int_F)[p]
        print(paste0("Covariates - Females: ", pheno))
        # Create a temporary dataframe: join all cov with the p phenotype
        temp <- left_join(pheno_int_F[, c(1,p)], cov_F, by = "IID")
	    # change name of p phenotype (col 2 in temp) to "pheno"
        colnames(temp)[2] <- "pheno"
        # Regress out covariates
        formula <- as.formula(paste0("pheno ~ age + age2 + batch + east_1km + north_1km + ", pcs_sum))
        pheno_int_cor_F[, pheno] <- residuals(lm(formula, data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_F, paste0(OUT_DIR, "/pheno_continuous_test_INT_age_age2_batch_location_PCs_F.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

