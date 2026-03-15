# Assortative mating null model (shuffling alpha, beta and gamma) for CNV-trait pairs

################################################################################
# Libraries
################################################################################

install.packages("data.table")
install.packages("dplyr")
library(data.table)
library(dplyr)

################################################################################
# Parameters
################################################################################

# Parse arguments if provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) coeff <- args[1]
if (length(args) >= 2) trials <- as.integer(args[2])
cat("coeff =", coeff, "\n")
cat("trials =", trials, "\n")

file_path_result <- paste0(coeff, ".shuffling.", trials,".txt")

################################################################################
# STEP 1: Load input files
################################################################################

# 1. Phenotype list
category_df <- read.table('/mnt/project/data/unique_phenotypes_annotated.43.txt', sep='\t', header=TRUE)
ordered_phenotypes_abbrev <- category_df[['PHENO']]

# 2. Phenotype, PGS_trans and CNV data
pheno_pgs_cnv <- fread('/mnt/project/data/pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.csv')

# 3. Couples
couples <- fread('/mnt/project/data/couples.test.csv')

# Add phenotype to couples 
pheno_F <- pheno_pgs_cnv %>% select(all_of(c("IID", ordered_phenotypes_abbrev)))
colnames(pheno_F) <- paste0(colnames(pheno_F), "_F")
pheno_M <- pheno_pgs_cnv %>% select(all_of(c("IID", ordered_phenotypes_abbrev)))
colnames(pheno_M) <- paste0(colnames(pheno_M), "_M")

pheno_couples <- merge(couples,       pheno_F, by = "IID_F",  all.x = TRUE)
pheno_couples <- merge(pheno_couples, pheno_M, by = "IID_M",  all.x = TRUE)

# 4. CNV-GWAS
# Import data on CNV signals (CNVR associated with a phenotype) and select only continuous phenotypes
CNVR <- read.table("/mnt/project/data/CNV_GWAS_signals_ranking_v1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Select only continuous phenotypes, not sex-specific + remove chromosome X
CNVR <- subset(CNVR, TYPE == "continuous" & SEX == "All" & CHR != 'X')
CNVR$CHR <- as.integer(CNVR$CHR)
CNVR$CNVR_START <- as.integer(CNVR$CNVR_START)
CNVR$CNVR_STOP <- as.integer(CNVR$CNVR_STOP)
CNVR$CNV_REGION <- paste(CNVR$CHR, CNVR$CNVR_START, CNVR$CNVR_STOP, sep = "_")

CNVR$PGS_START <- as.integer(CNVR$CNVR_START) - 250000
CNVR$PGS_STOP <- as.integer(CNVR$CNVR_STOP) + 250000
CNVR$PGS_REGION <- paste(CNVR$CHR, CNVR$PGS_START, CNVR$PGS_STOP, sep = "_")

# Order by phenotype
CNVR <- CNVR[order(CNVR$PHENO, CNVR$CHR, CNVR$CNVR_START), ]
# There is one top model which is both mirror (M) and duplication (D). Use M since it's more general.
both_model_i <- which(CNVR$TOP_MODEL == "M-DUP")
CNVR[both_model_i, "TOP_MODEL"] <- "M"
stopifnot(all(CNVR$TOP_MODEL != "M-DUP"))

# Some checks
stopifnot(nrow(CNVR) == 119)
stopifnot(length(unique(CNVR[['PHENO']])) == 43)

################################################################################
# STEP 2: Compute each correlation in the AM path
################################################################################

# Create empty data frames for saving correlations among couples and within individuals (1 row per CNVR)
corr_df <- data.frame(
    `pheno_F-pheno_M` = numeric(nrow(CNVR)),
    `pheno-PGS`       = numeric(nrow(CNVR)),
    `pheno-CNV`       = numeric(nrow(CNVR)),
    `PGS-CNV.total`   = numeric(nrow(CNVR)),
    row.names         = paste0(CNVR$PHENO, "_CNV_", CNVR$CNV_REGION),
    check.names = FALSE
)

# Iterate through each CNVR-trait pair
for (i in seq_len(nrow(CNVR))) {
    pheno <- CNVR[i, "PHENO"]   
    pheno_CNVR <- paste0(pheno, "_CNV_", CNVR[i, "CNV_REGION"])
    pheno_PGSR <- paste0(pheno, "_chr", CNVR[i, "PGS_REGION"], "_PGS")
    if(!pheno_CNVR %in% colnames(pheno_pgs_cnv)) stop("Missing CNV col: ", pheno_CNVR)
    if(!pheno_PGSR %in% colnames(pheno_pgs_cnv)) stop("Missing PGS col: ", pheno_PGSR)

    # Compute correlations
    # 'pheno_F-pheno_M'
    valid_data <- pheno_couples[!is.na(pheno_couples[[paste0(pheno, "_F")]]) & !is.na(pheno_couples[[paste0(pheno, "_M")]]), ]
    if (nrow(valid_data) > 1) {
        test <- cor.test(valid_data[[paste0(pheno, "_F")]], valid_data[[paste0(pheno, "_M")]], method = "pearson")
        corr_df[pheno_CNVR, "pheno_F-pheno_M"] <- unname(test$estimate)
    }

    # 'pheno-PGS_trans'
    valid_data <- pheno_pgs_cnv[!is.na(pheno_pgs_cnv[[pheno]]) & !is.na(pheno_pgs_cnv[[pheno_PGSR]]), ]
    if (nrow(valid_data) > 1) {
        test <- cor.test(valid_data[[pheno]], valid_data[[pheno_PGSR]], method = "pearson")
        corr_df[pheno_CNVR, "pheno-PGS"] <- unname(test$estimate)
    }

    # 'pheno-CNV'
    valid_data <- pheno_pgs_cnv[!is.na(pheno_pgs_cnv[[pheno]]) & !is.na(pheno_pgs_cnv[[pheno_CNVR]]), ]
    if (nrow(valid_data) > 1) {
        test <- cor.test(valid_data[[pheno]], valid_data[[pheno_CNVR]], method = "pearson")
        corr_df[pheno_CNVR, "pheno-CNV"] <- unname(test$estimate)
    }

    # 'pgs-CNV' ("total")
    valid_data <- pheno_pgs_cnv[!is.na(pheno_pgs_cnv[[pheno_PGSR]]) & !is.na(pheno_pgs_cnv[[pheno_CNVR]]), ]
    if (nrow(valid_data) > 1) {
        test <- cor.test(valid_data[[pheno_PGSR]], valid_data[[pheno_CNVR]], method = "pearson")
        corr_df[pheno_CNVR, "PGS-CNV.total"] <- unname(test$estimate)
    }
}


################################################################################
# STEP 3: Shuffle alpha, beta or gamma
################################################################################

# Get lists of CNV and PGS names
set.seed(81)

if (coeff=="alpha") {
    cnv_list <- c()
    for (i in seq_len(nrow(CNVR))) {
        pheno <- CNVR[i, "PHENO"]
        pheno_CNVR <- paste0(pheno, "_CNV_", CNVR[i, "CNV_REGION"])
        cnv_list <- c(cnv_list, pheno_CNVR)
    }

}

if (coeff=="gamma") {
    pgs_list <- c()
    for (i in seq_len(nrow(CNVR))) {
        pheno <- CNVR[i, "PHENO"]
        pheno_PGSR <- paste0(pheno, "_chr", CNVR[i, "PGS_REGION"], "_PGS")
        pgs_list <- c(pgs_list, pheno_PGSR)
    }
}

# Initialize vector to store indirect-total CNV-PGS correlations for randomized coefficients
ind_tot_corr <- numeric(trials)

for (trial in seq_len(trials)){
    # Do the shuffling 
    if (coeff=="alpha") {
        cnv_list <- sample(cnv_list)

    } else if (coeff=="gamma") {
        pgs_list <- sample(pgs_list)

    } else if (coeff=="beta") {
        couples_shuffled <- couples
        couples_shuffled$IID_M <- sample(couples_shuffled$IID_M) # shuffle the male partner
        pheno_couples_shuffled <- merge(couples_shuffled,       pheno_F, by = "IID_F",  all.x = TRUE)
        pheno_couples_shuffled <- merge(pheno_couples_shuffled, pheno_M, by = "IID_M",  all.x = TRUE)

    } else {
        print('Option not valid')
    }

    for (j in seq_len(nrow(CNVR))) {
        pheno <- CNVR[j, "PHENO"]
        pheno_CNVR <- paste0(pheno, "_CNV_", CNVR[j, "CNV_REGION"])
        pheno_PGSR <- paste0(pheno, "_chr", CNVR[j, "PGS_REGION"], "_PGS")
        if(!pheno_CNVR %in% colnames(pheno_pgs_cnv)) stop("Missing CNV col: ", pheno_CNVR)
        if(!pheno_PGSR %in% colnames(pheno_pgs_cnv)) stop("Missing PGS col: ", pheno_PGSR)

        if (coeff=="alpha"){
            # 'pheno-CNV'
            # Assign a random CNVR (j)
            valid_data <- pheno_pgs_cnv[!is.na(pheno_pgs_cnv[[pheno]]) & !is.na(pheno_pgs_cnv[[cnv_list[j]]]), ]
            if (nrow(valid_data) > 1) {
                test <- cor.test(valid_data[[pheno]], valid_data[[cnv_list[j]]], method = "pearson")
                corr_df[pheno_CNVR, "pheno-CNV"] <- unname(test$estimate)
            }
        } else if (coeff=="gamma"){
            # 'pheno-PGS_trans'
            # Assign a random PGS (j)
            valid_data <- pheno_pgs_cnv[!is.na(pheno_pgs_cnv[[pheno]]) & !is.na(pheno_pgs_cnv[[pgs_list[j]]]), ]
            if (nrow(valid_data) > 1) {
                test <- cor.test(valid_data[[pheno]], valid_data[[pgs_list[j]]], method = "pearson")
                corr_df[pheno_CNVR, "pheno-PGS"] <- unname(test$estimate)
            }
        } else if (coeff=="beta") {
            # 'pheno_F-pheno_M'
            # Use random couples
            valid_data <- pheno_couples_shuffled[!is.na(pheno_couples_shuffled[[paste0(pheno, "_F")]]) & !is.na(pheno_couples_shuffled[[paste0(pheno, "_M")]]), ]
            if (nrow(valid_data) > 1) {
                test <- cor.test(valid_data[[paste0(pheno, "_F")]], valid_data[[paste0(pheno, "_M")]], method = "pearson")
                corr_df[pheno_CNVR, "pheno_F-pheno_M"] <- unname(test$estimate)
            }
        }
    }

    # Compute new indirect PGS-CNV correlation
    corr_df[['PGS-CNV.indirect']] <- corr_df[['pheno-CNV']]*corr_df[['pheno_F-pheno_M']]*corr_df[['pheno-PGS']]
    # Compute new indirect vs total PGS-CNV correlation
    ind_tot_corr[trial] <- unname(cor.test(corr_df[["PGS-CNV.total"]], corr_df[["PGS-CNV.indirect"]])$estimate)
    print(paste(trial, ind_tot_corr[trial]))
}

# Save corr(indirect CNV-PGS corr, direct CNV-PGS corr)
ind_tot_corr_df <- data.frame(ind_tot_corr)
write.table(ind_tot_corr_df, file = file_path_result, row.names = FALSE, col.names = FALSE, quote = FALSE)

