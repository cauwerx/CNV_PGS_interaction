library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(data.table)

#------------------------------------------------------------------------------#
# Parameters
#------------------------------------------------------------------------------#

args = commandArgs(trailingOnly=T)
WORK="/data/LDpred2/"
SCRATCH="/scratch/LDpred2_train/"

# input
pheno=args[1]
arg_p=as.numeric(args[2])
multi_auto_rds_p = args[3]
sum_s_out = args[4]
map_p = paste0(WORK, "Data/map_hm3_plus.rds") # use hapmap+
sum_stat_p = paste0("/data/UKBB/GWAS/Regenie/step2/ukb.train_", pheno, ".regenie.gz")

# output
p_value = seq_log(1e-4, 0.2, length.out = 10)[arg_p] 
tmp_p = paste0(SCRATCH, "tmp-data-", pheno, ".", arg_p)
cat("Proportion used:", p_value, '\n')

#------------------------------------------------------------------------------#
# STEP 1: Obtain HapMap3 SNPs
#------------------------------------------------------------------------------#
info = readRDS(map_p)
cat("HapMap3 SNPs imported\n")

#------------------------------------------------------------------------------#
# STEP 2: Load the summary statistic file
#------------------------------------------------------------------------------#
maf_filter=0.01; info_filter=0.8

sumstats <- bigreadr::fread2(sum_stat_p, select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE"))
names(sumstats) = c("chr", "pos", "rsid", "a0", "a1", "a1f", "INFO", "n_eff", "beta", "beta_se")
# Use civilised rsid (chr:pos_ref_alt) because there are duplicated variants
sumstats$civ_rsid = paste0(sumstats$chr, ":", sumstats$pos, "_", sumstats$a1, "_", sumstats$a0)
sumstats$rev_civ_rsid = paste0(sumstats$chr, ":", sumstats$pos, "_", sumstats$a0, "_", sumstats$a1)
info$civ_rsid = paste0(info$chr, ":", info$pos, "_", info$a1, "_", info$a0)
# Remove ambigous SNPs
ambigous_snps = ( (sumstats$a0=="A" & sumstats$a1=="T") | (sumstats$a0=="T" & sumstats$a1=="A") | (sumstats$a0=="G" & sumstats$a1=="C") | (sumstats$a0=="C" & sumstats$a1=="G"))
cat("#ambigous snps:", sum(ambigous_snps), "\n")
cat("#non ambigous snps:", sum(!ambigous_snps), "\n")
sumstats_filtered = sumstats[! ambigous_snps, ]
# Check that no rsid are duplicated
stopifnot(length(unique(sumstats_filtered$civ_rsid)) == nrow(sumstats_filtered))
# Filter out sumstats_snps
sumstats_filtered <- sumstats_filtered[sumstats_filtered$civ_rsid%in% info$civ_rsid | sumstats_filtered$rev_civ_rsid %in% info$civ_rsid,]
sumstats_filtered_maf_info = sumstats_filtered[sumstats_filtered$a1f > maf_filter & sumstats_filtered$a1f < 1-maf_filter & sumstats_filtered$INFO > info_filter,]
cat("Summary statistic imported and filtered.\n")

# Validation

#1. check if any ambigous snps are left
ambigous_snps_filtered_maf_info = ( (sumstats_filtered_maf_info$a0=="A" & sumstats_filtered_maf_info$a1=="T") | (sumstats_filtered_maf_info$a0=="T" & sumstats_filtered_maf_info$a1=="A") | (sumstats_filtered_maf_info$a0=="G" & sumstats_filtered_maf_info$a1=="C") | (sumstats_filtered_maf_info$a0=="C" & sumstats_filtered_maf_info$a1=="G"))
stopifnot(sum(ambigous_snps_filtered_maf_info)==0)
#2. check that there are no info < 0.8
stopifnot(all(sumstats_filtered_maf_info$INFO > info_filter))
#3. check that there are no MAF < 0.01
stopifnot(all(sumstats_filtered_maf_info$a1f < info_filter | sumstats_filtered_maf_info$a1f > 1-info_filter))
#4. GWAS QC
sd_Gj_est  <- sqrt(sumstats_filtered_maf_info$n_eff * sumstats_filtered_maf_info$beta_se^2 + sumstats_filtered_maf_info$beta^2)
sd_y <- quantile(sqrt(0.5 * (sumstats_filtered_maf_info$n_eff * sumstats_filtered_maf_info$beta_se^2 + sumstats_filtered_maf_info$beta^2)), probs = 0.01, na.rm = TRUE)
sd_s <- sd_y / sd_Gj_est
sd_ldref <- sqrt(2 * sumstats_filtered_maf_info$a1f * (1 - sumstats_filtered_maf_info$a1f) * sumstats_filtered_maf_info$INFO)
is_bad <- sd_s < (0.5 * sd_ldref) | sd_s > (sd_ldref + 0.1) | sd_s < 0.05 | sd_ldref < 0.05
stopifnot(all(!is_bad))

#------------------------------------------------------------------------------#
# STEP 3: Get pre-computed LD
#------------------------------------------------------------------------------#

# Initialize variables
ld <- NULL
corr <- NULL

# Loop through chromosomes
tmp <- tempfile(tmpdir = tmp_p)
for (chr in 1:22) {
	cat(paste("Processing chromosome", chr, "\n"))

	# Load the precomputed LD matrix for this chromosome
	ld_p = paste0(WORK, "Data/ldref_plus/LD_with_blocks_chr", chr,".rds") # use the hapmap + LD
	corr0 <- readRDS(ld_p)  # Load precomputed LD matrix

	# Combine into a single SFBM object
	if (chr == 1) {
	    ld <- Matrix::colSums(corr0^2)  # LD scores for first chromosome
	    corr <- as_SFBM(corr0, compact = TRUE)  # Initialize the SFBM
	} else {
	    ld <- c(ld, Matrix::colSums(corr0^2))  # Append LD scores
	    corr$add_columns(corr0, nrow(corr))  # Append columns to the SFBM
	}
}

stopifnot(length(ld) == nrow(info))
stopifnot(length(ld) == ncol(corr))
stopifnot(length(ld) == nrow(corr))
cat("Pre-computed LD imported and filtered.\n")

#------------------------------------------------------------------------------#
# STEP 4: Run LDpred2
#------------------------------------------------------------------------------#

# Estimate of h2 from LD Score regression
stopifnot(all(sumstats_filtered_maf_info$civ_rsid %in% info$civ_rsid | sumstats_filtered_maf_info$rev_civ_rsid %in% info$civ_rsid)) # do we only have matching positions in sumstats_filtered_maf_info?

matched_variants = info$rsid[info$civ_rsid %in% sumstats_filtered_maf_info$civ_rsid  | info$civ_rsid %in% sumstats_filtered_maf_info$rev_civ_rsid]
stopifnot(length(matched_variants) == nrow(sumstats_filtered_maf_info))
ind.ref = match(matched_variants, info$rsid) 
ind.df_beta = match(matched_variants, sumstats_filtered_maf_info$rsid) 

stopifnot(all(sumstats_filtered_maf_info[ind.df_beta, ]$rsid == info[ind.ref, ]$rsid ))
ld_s_new <- ld[ind.ref]
sum_s_new = sumstats_filtered_maf_info[ind.df_beta,]

# save positions of sum_s_new
output_table = data.frame(
	"rsid" = sum_s_new$rsid,
	"chr_name" = sum_s_new$chr,
	"chr_position" = sum_s_new$pos,
	"effect_allele" = sum_s_new$a1,
	"other_allele" = sum_s_new$a0
	)

# only write output table one time
if(arg_p==1) {write.table(output_table, file=sum_s_out, quote=F, row.names=F, sep='\t')}

ldsc <- with(sum_s_new, snp_ldsc(ld_s_new, nrow(info), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))
ldsc_h2_est <- ldsc[["h2"]]
cat("H^2 estimated:", ldsc_h2_est, '\n')

(NCORES <- nb_cores())
cat(paste("Ncores =", NCORES, "\n"))
coef_shrink = 0.95

multi_auto <- snp_ldpred2_auto(
	corr, sum_s_new, h2_init = ldsc_h2_est, ind.corr = ind.ref,
	vec_p_init = c(p_value), ncores = NCORES,
	burn_in = 500, num_iter = 500, report_step = 20,
	allow_jump_sign = FALSE, use_MLE = TRUE, shrink_corr = coef_shrink)

saveRDS(multi_auto, multi_auto_rds_p, compress=T)
cat(paste("Output written in: ", multi_auto_rds_p, "\n"))

