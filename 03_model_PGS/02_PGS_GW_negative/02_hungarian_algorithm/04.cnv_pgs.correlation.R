# Computes CNV-PGS correlation for both empirical and null CNV-PGS links

library(data.table)
setwd('/penetrance_at_scale/data/')

################################################################################
# STEP 1: Load data
################################################################################

# 1. Load empirical and null sets
assign_table <- read.csv("07_sensitivity_analyses/02.trait_cnv.empirical_and_null_sets.csv", header = TRUE, stringsAsFactors = FALSE)

# 2. Load genotype and phenotype data
geno_pheno <- as.data.frame(fread("06_pheno_geno_summary/pheno_pgs_cnv.autosomes.non_sex_specific.csv"))

stopifnot(all(assign_table$PGS %in% names(geno_pheno)))
stopifnot(all(assign_table$CNV %in% names(geno_pheno)))
stopifnot(all(assign_table$new_CNV %in% names(geno_pheno)))

################################################################################
# STEP 2: Compute CNV-PGS correlations
################################################################################

assign_table$PGS_CNV.r <- mapply(
    FUN = function(pgs, cnv) cor(geno_pheno[[pgs]], geno_pheno[[cnv]], use = "complete.obs"),
    assign_table$PGS,
    assign_table$CNV
)

assign_table$PGS_newCNV.r <- mapply(
    FUN = function(pgs, cnv) cor(geno_pheno[[pgs]], geno_pheno[[cnv]], use = "complete.obs"),
    assign_table$PGS,
    assign_table$new_CNV
)

assign_table$PGS_CNV.r_sign <- assign_table$PGS_CNV.r * sign(assign_table$trait_CNV.r)
assign_table$PGS_newCNV.r_sign <- assign_table$PGS_newCNV.r * sign(assign_table$trait_newCNV.r)

cat('Mean of empirical CNV-PGS correlations with sign-adjustment:', mean(assign_table$PGS_CNV.r_sign), '\n')
cat('Mean of null CNV-PGS correlations with sign-sdjustment:', mean(assign_table$PGS_newCNV.r_sign), '\n')

################################################################################
# STEP 3: Perform two-sided paired t-test to compare the groups
################################################################################

cat("t-test p-value:", t.test(assign_table$PGS_CNV.r_sign, assign_table$PGS_newCNV.r_sign, paired=TRUE)$p.value, '\n')

################################################################################
# STEP 4: Save results
################################################################################

fwrite(assign_table, "07_sensitivity_analyses/03.pgs_cnv.corr.csv")

