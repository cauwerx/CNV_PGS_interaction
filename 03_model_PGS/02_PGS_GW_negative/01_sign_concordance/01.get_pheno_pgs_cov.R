# Saves parameters from real data for sign-concordance simulations
# Specifically: PGS and pheno covariance matrices,
# and effect sizes for the 2 models: pheno~CNV and PGS~CNV

#############################################################################################
# Libraries
#############################################################################################
library(data.table)
library(dplyr)
project_DIR = "/mnt/project/"

#############################################################################################
# Load data
#############################################################################################

# Phenotype and genotype information
pheno_pgs_cnvs <- as.data.frame(
    fread(file = file.path(project_DIR, '/data/pheno_pgs_cnv.autosomes.non_sex_specific.csv'))
)
# Load model results results
CNV_freq <- fread(file = file.path(project_DIR, '/sign_concordance/data/CNV_freq.txt'))
additive_model <- fread(file = file.path(project_DIR, '/sign_concordance/data/full_model.txt'))
additive_model <- merge(additive_model, CNV_freq, by = c("ID", "phenotype", "cytogenic_band", "top_model"))
# Assign most general model if more than one is significant
additive_model[additive_model$top_model=='M-DUP']$top_model <- 'M'

#############################################################################################
# STEP 1: Perform both regressions and save variant info
#############################################################################################

param_df <- additive_model[,c('ID', 'phenotype', 'cytogenic_band', 'chromosome', 'start_CNVR', 'end_CNVR', 'top_model')]
param_df$n <- integer(0)
param_df$CNV_n <- integer(0)
param_df$pheno_sd <- numeric(0)
param_df$PGS_mean <- numeric(0)
param_df$PGS_sd <- numeric(0)
param_df$CNV_coefficient <- numeric(0)
param_df$PGS_coefficient <- numeric(0)

for (i in 1:nrow(additive_model)){
    # Get model
    pheno_name <- additive_model[[i, 'phenotype']]
    pgs_name <- paste0(pheno_name, '_PGS')
    cnv_name <- paste0(pheno_name, '_CNV_', additive_model[i, 'chromosome'], '_', additive_model[i, 'start_CNVR'], '_', additive_model[i, 'end_CNVR'])
    top_model <- additive_model[i, 'top_model']

    formula <- as.formula(paste(pheno_name, '~', pgs_name, '+', cnv_name))
    model <- lm(formula, data=pheno_pgs_cnvs)

    # Check I get same coefficients
    # PGS effect
    PGS_coeff <- model$coefficients[pgs_name]
    stopifnot(round(PGS_coeff, 9) == round(additive_model[[i, 'effect_PGS_GW']], 9))
    # CNV effect
    CNV_coeff <- unname(model$coefficients[cnv_name])
    stopifnot(round(CNV_coeff, 9) == round(additive_model[[i, 'effect_CNV']], 9))

    # Get n (individuals with non NA for pheno, PGS and CNV)
    n <- nrow(model$model)

    # Get phenotype_sd
    pheno_sd <- sd(model$model[[pheno_name]])

    # Get mean and SD for PGS
    PGS_mean <- mean(model$model[[pgs_name]])
    PGS_sd <- sd(model$model[[pgs_name]])

    # Get CNV n (NA removed)
    CNV_n <- sum(model$model[[cnv_name]] != 0)

    # Get CNV n old, it may include individuals with missing pheno info
    if (top_model=="M"){
        CNV_n_old <- additive_model[i, 'number_CNV']
    } else if (top_model=="DUP") {
        CNV_n_old <- additive_model[i, 'number_duplication']
    } else if (top_model=="DEL") {
        CNV_n_old <- additive_model[i, 'number_deletion']
    } else {
        print("Model not available.")
    }

    stopifnot(CNV_n <= CNV_n_old)
    print(paste(cnv_name, CNV_n, CNV_n_old))

    # Save results
    param_df[i, 'n'] <- n
    param_df[i, 'CNV_n'] <- CNV_n
    param_df[i, 'pheno_sd'] <- pheno_sd
    param_df[i, 'PGS_mean'] <- PGS_mean
    param_df[i, 'PGS_sd'] <- PGS_sd
    param_df[i, 'CNV_coefficient'] <- CNV_coeff
    param_df[i, 'PGS_coefficient'] <- PGS_coeff
}

fwrite(param_df, 'additive_model.parameters.csv')

#############################################################################################
# STEP 2: Get phenotype and PGS covariance matrices
#############################################################################################

unique_pheno <- unique(additive_model[['phenotype']])

pheno_df <- pheno_pgs_cnvs[unique_pheno]
pheno_cov <- data.frame(cov(pheno_df, use="pairwise.complete.obs"))

pgs_df <- pheno_pgs_cnvs[paste0(unique_pheno, "_PGS")]
pgs_cov <- data.frame(cov(pgs_df, use="pairwise.complete.obs"))

fwrite(pheno_cov, 'pheno.cov_matrix.csv', row.names=TRUE)
fwrite(pgs_cov, 'pgs.cov_matrix.csv', row.names=TRUE)


