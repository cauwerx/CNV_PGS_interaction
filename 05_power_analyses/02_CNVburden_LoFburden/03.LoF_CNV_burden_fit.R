# Fit normal, negative binomial and Poisson models to CNV_burden and LoF_burden

###########################################################################################
# Libraries
###########################################################################################
library(data.table)
library(MASS)
library(ggplot2)

###########################################################################################
# Load data
###########################################################################################

# Get white-british individuals
sample_qc <- as.data.frame(fread("/mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb_sqc_v2.txt", header = F, select = c(3:68)))
colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")
sample_eid <- as.data.frame(fread("/mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c1_b0_v2.fam", header = F, select = c(1,2,5), col.names = c("fid", "eid", "sex")))
df <- cbind(sample_eid, sample_qc)
WB_df <- df[df$white_british == 1, "eid", drop = FALSE]

# Load CNV burden and LoF burden
cnv_burden = as.data.frame(fread('/mnt/project/data/CNV_burden/CNV_burden.txt.gz', sep='\t'))
lof_burden = as.data.frame(fread('/mnt/project/data/LoF_burden/LoF_burden.total.csv', sep=','))

# Filter by WB individuals
cnv_burden = merge(WB_df, cnv_burden, by.x = "eid", by.y = "IID", all.x = TRUE)
lof_burden = merge(WB_df, lof_burden, by.x = "eid", by.y = "IID", all.x = TRUE)
burden = merge(lof_burden, cnv_burden, by='eid', all = TRUE)

###########################################################################################
# STEP 1: Compute LoF burden - CNV burden correlation
###########################################################################################
mask <- !is.na(burden$LoF_total) & !is.na(burden$BURDEN_GENES)
valid_data <- burden[mask,]
WB_n <- nrow(valid_data)
print(cor.test(valid_data$LoF_total, valid_data$BURDEN_GENES))

###########################################################################################
# STEP 2: Fit models to CNV burden and LoF burden data
###########################################################################################
set.seed(123)

# Get parameters

# -- CNV burden --
# Normal
norm_fit_cnv  <- fitdistr(burden[['BURDEN_GENES']][!is.na(burden[['BURDEN_GENES']])], "normal")
norm_sd_cnv   <- unname(norm_fit_cnv$estimate["sd"])
norm_mean_cnv <- unname(norm_fit_cnv$estimate["mean"])
# Negative Binomial
nb_fit_cnv  <- fitdistr(burden[['BURDEN_GENES']][!is.na(burden[['BURDEN_GENES']])], "Negative Binomial")
nb_size_cnv <- unname(nb_fit_cnv$estimate["size"])
nb_mu_cnv   <- unname(nb_fit_cnv$estimate["mu"])
# Poisson
pois_fit_cnv    <- fitdistr(burden[['BURDEN_GENES']][!is.na(burden[['BURDEN_GENES']])], "Poisson")
pois_lambda_cnv <- unname(pois_fit_cnv$estimate["lambda"])

# -- LoF burden --
# Normal
norm_fit_lof  <- fitdistr(burden[['LoF_total']][!is.na(burden[['LoF_total']])], "normal")
norm_sd_lof   <- unname(norm_fit_lof$estimate["sd"])
norm_mean_lof <- unname(norm_fit_lof$estimate["mean"])
# Negative Binomial
nb_fit_lof  <- fitdistr(burden[['LoF_total']][!is.na(burden[['LoF_total']])], "Negative Binomial")
nb_size_lof <- unname(nb_fit_lof$estimate["size"])
nb_mu_lof   <- unname(nb_fit_lof$estimate["mu"])
# Poisson
pois_fit_lof    <- fitdistr(burden[['LoF_total']][!is.na(burden[['LoF_total']])], "Poisson")
pois_lambda_lof <- unname(pois_fit_lof$estimate["lambda"])


###########################################################################################
# STEP 3: Simulate CNV burden and LoF burden
###########################################################################################

normal_cnv <- rnorm(n=WB_n, mean=norm_mean_cnv, sd=norm_sd_cnv)
normal_lof <- rnorm(n=WB_n, mean=norm_mean_lof, sd=norm_sd_lof)
nbinom_cnv <- rnbinom(n=WB_n, size=nb_size_cnv, mu=nb_mu_cnv)
nbinom_lof <- rnbinom(n=WB_n, size=nb_size_lof, mu=nb_mu_lof)
poisson_cnv <- rpois(n=WB_n, lambda=pois_lambda_cnv)
poisson_lof <- rpois(n=WB_n, lambda=pois_lambda_lof)

valid_data$normal_cnv <- normal_cnv
valid_data$normal_lof <- normal_lof
valid_data$nbinom_cnv <- nbinom_cnv
valid_data$nbinom_lof <- nbinom_lof
valid_data$poisson_cnv <- poisson_cnv
valid_data$poisson_lof <- poisson_lof

dist_cnv_df <- valid_data[, c('BURDEN_GENES', 'normal_cnv', 'nbinom_cnv', 'poisson_cnv')]
dist_lof_df <- valid_data[, c('LoF_total', 'normal_lof', 'nbinom_lof', 'poisson_lof')]

dist_cnv_df[] <- sapply(dist_cnv_df, as.numeric)
dist_lof_df[] <- sapply(dist_lof_df, as.numeric)

dist_cnv_longdf <- melt(setDT(dist_cnv_df),
                        measure.vars=c('BURDEN_GENES', 'normal_cnv', 'nbinom_cnv', 'poisson_cnv'),
                        variable.name='type',
                        value.name="value")

dist_lof_longdf <- melt(setDT(dist_lof_df),
                        measure.vars=c('LoF_total', 'normal_lof', 'nbinom_lof', 'poisson_lof'),
                        variable.name='type',
                        value.name="value")

