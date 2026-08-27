# Get the Pearson correlation coefficients for all trait-CNV_encoding combinations
library(data.table)

setwd('/penetrance_at_scale/data/')

################################################################################
# STEP 1: Load data
################################################################################

# 1. Load CNV-GWAS signals
cnv_gwas <- read.table("./CNV_GWAS_signals_ranking_v1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Filter for continuous phenotypes, only autosomal CHR and sex=All
cnv_gwas <- subset(cnv_gwas, TYPE == "continuous" & CHR != 'X' & SEX=='All')

# There is one row with two top models, assign the most general
cnv_gwas$TOP_MODEL[cnv_gwas$TOP_MODEL == "M-DUP"] <- "M"

# CNV names
cnv_gwas$CNV <- paste(cnv_gwas$PHENO, 'CNV', cnv_gwas$CHR, cnv_gwas$CNVR_START, cnv_gwas$CNVR_STOP, sep='_')
cnv_gwas[c("CHR", "CNVR_START", "CNVR_STOP")] <- lapply(cnv_gwas[c("CHR", "CNVR_START", "CNVR_STOP")], FUN=as.integer)
cnv_gwas <- cnv_gwas[order(cnv_gwas$CHR, cnv_gwas$CNVR_START),]
n <- nrow(cnv_gwas)
stopifnot(n==119)
unique_traits <- sort(unique(cnv_gwas$PHENO))
unique_cnvs <- unique(cnv_gwas$CNV)
trait_n <- length(unique_traits)
cnv_n <- length(unique_cnvs)
stopifnot(n == cnv_n)

# 2. Load genotype and phenotype data
geno_pheno <- as.data.frame(fread("06_pheno_geno_summary/pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.csv"))

stopifnot(all(unique_traits %in% names(geno_pheno)))
stopifnot(all(unique_cnvs %in% names(geno_pheno)))

################################################################################
# STEP 2: Create trait-CNV matrix of Pearson correlation coefficients
################################################################################

r_m <- matrix(NA_real_, nrow=trait_n, ncol=cnv_n)
rownames(r_m) <- unique_traits
colnames(r_m) <- unique_cnvs

for (trait in unique_traits){
  for (cnv in unique_cnvs){
    r <- cor.test(geno_pheno[[trait]], geno_pheno[[cnv]], use = "complete.obs")$estimate
    r_m[trait, cnv] <- r
    cat(trait, cnv, r, '\n')

  }
}

fwrite(as.data.table(r_m, keep.rownames = "trait"), "07_sensitivity_analyses/01.trait_cnv.corr_coef.csv")


