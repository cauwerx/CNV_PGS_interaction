# Assortative mating null model (shuffling alpha, beta and gamma) for CNV burden

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

file_path_result <- paste0(coeff, ".burden.shuffling.", trials,".txt")

################################################################################
# STEP 1: Load input files
################################################################################

# 1. Phenotype list
category_df <- read.table('/mnt/project/data/unique_phenotypes_annotated.43.txt', sep='\t', header=TRUE)
ordered_phenotypes_abbrev <- category_df[['PHENO']]

# 2. Phenotype, PGS_GW and CNV data
pheno_pgs_cnv <- as.data.frame(fread('/mnt/project/data/pheno_pgs_cnv.autosomes.non_sex_specific.csv'))

# 3. Couples
couples <- as.data.frame(fread('/mnt/project/data/couples.test.csv'))

# Add phenotype to couples df
pheno_F <- pheno_pgs_cnv %>% select(all_of(c("IID", ordered_phenotypes_abbrev)))
colnames(pheno_F) <- paste0(colnames(pheno_F), "_F")
pheno_M <- pheno_pgs_cnv %>% select(all_of(c("IID", ordered_phenotypes_abbrev)))
colnames(pheno_M) <- paste0(colnames(pheno_M), "_M")

pheno_couples <- merge(couples,       pheno_F, by = "IID_F",  all.x = TRUE)
pheno_couples <- merge(pheno_couples, pheno_M, by = "IID_M",  all.x = TRUE)

# 4. CNV_burden
cnv_burden <- as.data.frame(fread('/mnt/project/data/CNV_burden/CNV_burden.txt.gz'))

# Remove individual CNVs
pheno_pgs_cnvburden <- pheno_pgs_cnv[, !grepl("_CNV_", names(pheno_pgs_cnv))]
# Add CNV_burden
pheno_pgs_cnvburden <- merge(pheno_pgs_cnvburden, cnv_burden[,c('IID', 'BURDEN_GENES')], by='IID')

# 5. pheno-CNV_burden associations
pheno_affected <- c("neuroticism", "GS", "vitamin_D", "WHR", "fluid_intelligence", "WHRadjBMI", "heart_rate", "neutrophil_count", "FVC",  "body_fat_mass", "WBC_count", "BMI", "CRP", "ApoA", "HbA1c", "cystatinC", "BMD", "HDL", "ALP", "platelet_count", "height")


################################################################################
# STEP 2: Compute each correlation in the AM path
################################################################################

# Create empty data frames for saving correlations among couples and within individuals (1 row per pheno)

pheno_num <- length(pheno_affected)

corr_df <- data.frame(
    `pheno_F-pheno_M`  = numeric(pheno_num),
    `pheno-PGS`        = numeric(pheno_num),
    `pheno-CNVburden`  = numeric(pheno_num),
    `PGS-CNVburden.total`    = numeric(pheno_num),
    `PGS-CNVburden.indirect` = numeric(pheno_num),
    row.names   = pheno_affected,
    check.names = FALSE
)

# Iterate through each pheno
for (pheno in pheno_affected) {
    pgs <- paste0(pheno, "_PGS")
    # Compute correlations
    # 'pheno_F-pheno_M'
    valid_data <- pheno_couples[!is.na(pheno_couples[[paste0(pheno, "_F")]]) & !is.na(pheno_couples[[paste0(pheno, "_M")]]), ]
    if (nrow(valid_data) > 1) {
        corr_df[pheno, "pheno_F-pheno_M"] <- cor(valid_data[[paste0(pheno, "_F")]], valid_data[[paste0(pheno, "_M")]])
    }

    # 'pheno-PGS_GW'
    valid_data <- pheno_pgs_cnvburden[!is.na(pheno_pgs_cnvburden[[pheno]]) & !is.na(pheno_pgs_cnvburden[[pgs]]), ]
    if (nrow(valid_data) > 1) {
        corr_df[pheno, "pheno-PGS"] <- cor(valid_data[[pheno]], valid_data[[pgs]])
    }

    # 'pheno-CNVburden'
    valid_data <- pheno_pgs_cnvburden[!is.na(pheno_pgs_cnvburden[[pheno]]) & !is.na(pheno_pgs_cnvburden[['BURDEN_GENES']]), ]
    if (nrow(valid_data) > 1) {
        corr_df[pheno, "pheno-CNVburden"] <- cor(valid_data[[pheno]], valid_data[['BURDEN_GENES']])
    }

    # 'PGS_GW-CNVburden' ("total")
    valid_data <- pheno_pgs_cnvburden[!is.na(pheno_pgs_cnvburden[[pgs]]) & !is.na(pheno_pgs_cnvburden[['BURDEN_GENES']]), ]
    if (nrow(valid_data) > 1) {
        corr_df[pheno, "PGS-CNVburden.total"] <- cor(valid_data[[pgs]], valid_data[['BURDEN_GENES']])
    }
}

################################################################################
# STEP 3: Shuffle alpha, beta or gamma
################################################################################

set.seed(81)

ind_tot_corr <- numeric(trials)

for (trial in seq_len(trials)){
    # Do the shuffling
    if (coeff=="alpha" | coeff=="gamma") {
        pheno_shuff_lst <- sample(pheno_affected) # shuffle pheno
    
    } else if (coeff=="beta") {
        couples_shuffled <- couples
        couples_shuffled$IID_M <- sample(couples_shuffled$IID_M) # shuffle the male partner
        pheno_couples_shuffled <- merge(couples_shuffled,       pheno_F, by = "IID_F",  all.x = TRUE)
        pheno_couples_shuffled <- merge(pheno_couples_shuffled, pheno_M, by = "IID_M",  all.x = TRUE)
        
    } else {
        print('Option not valid')
    }
    
    for (j in 1:length(pheno_affected)) {  
        pheno <- pheno_affected[j]
        pgs <- paste0(pheno, "_PGS")
        
        if (coeff=="alpha"){
            # 'pheno-CNVburden'
            # Assign a random pheno
            valid_data <- pheno_pgs_cnvburden[!is.na(pheno_pgs_cnvburden[[pheno_shuff_lst[j]]]) & !is.na(pheno_pgs_cnvburden[['BURDEN_GENES']]), ]
            if (nrow(valid_data) > 1) {
                corr_df[pheno, "pheno-CNVburden"] <- cor(valid_data[[pheno_shuff_lst[j]]], valid_data[['BURDEN_GENES']])
            }
        } else if (coeff=="gamma"){
            # 'pheno-PGS_trans'
            # Assign a random pheno
            valid_data <- pheno_pgs_cnvburden[!is.na(pheno_pgs_cnvburden[[pheno_shuff_lst[j]]]) & !is.na(pheno_pgs_cnvburden[[pgs]]), ]
            if (nrow(valid_data) > 1) {
                corr_df[pheno, "pheno-PGS"] <- cor(valid_data[[pheno_shuff_lst[j]]], valid_data[[pgs]])
            }
        } else if (coeff=="beta") {
            # 'pheno_F-pheno_M'
            # Use random couples
            valid_data <- pheno_couples_shuffled[!is.na(pheno_couples_shuffled[[paste0(pheno, "_F")]]) & !is.na(pheno_couples_shuffled[[paste0(pheno, "_M")]]), ]
            if (nrow(valid_data) > 1) {
                corr_df[pheno, "pheno_F-pheno_M"] <- cor(valid_data[[paste0(pheno, "_F")]], valid_data[[paste0(pheno, "_M")]])
            }    
        }
    }

    # Compute new indirect PGS-CNVburden correlation
    corr_df[['PGS-CNVburden.indirect']] <- corr_df[['pheno-CNVburden']]*corr_df[['pheno_F-pheno_M']]*corr_df[['pheno-PGS']]
    
    # Compute new indirect vs total PGS-CNVburden correlation
    ind_tot_corr[trial] <- cor(corr_df[["PGS-CNVburden.total"]], corr_df[["PGS-CNVburden.indirect"]])
}

# Save corr(indirect CNV-PGS corr, direct CNV-PGS corr)
ind_tot_corr_df <- data.frame(ind_tot_corr)
write.table(ind_tot_corr_df, file = file_path_result, row.names = FALSE, col.names = FALSE, quote = FALSE)


