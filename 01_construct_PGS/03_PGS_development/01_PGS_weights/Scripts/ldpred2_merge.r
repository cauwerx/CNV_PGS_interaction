library(bigsnpr)

# Merge the ldpred2 results generated for different values of p (the proportion of causal variants)

#------------------------------------------------------------------------------#
# Parameters
#------------------------------------------------------------------------------#
args = commandArgs(trailingOnly=T)
pheno = args[1]
ldpred_f = args[2] 
sum_s_p = args[3]
out_alpha = args[4]
out_beta = args[5]

WORK="/data/LDpred2/"
SCRATCH="/scratch/LDpred2_2025_05_02/"

# input
map_p = paste0(WORK, "Data/map_hm3_plus.rds") # use hapmap+

# output
cat(out_alpha, '\n')
cat(out_beta, '\n')

#------------------------------------------------------------------------------#
# STEP 1: Import the data
#------------------------------------------------------------------------------#
#rds_files = list.files(path=ldpred_f, pattern=paste0("^", pheno, ".", "*rds$"))
rds_files = paste0(ldpred_f, pheno, ".", c(1:10), ".multi_auto.rds")
print(rds_files)

multi_auto = list()
cat("Importing each proportion value.\n")
for(rds_idx in c(1:length(rds_files))){
	rds_p = rds_files[rds_idx]
	ldpred2_res = readRDS(rds_p)
	multi_auto[rds_idx] = ldpred2_res
	cat(rds_p, '\n')
}

#------------------------------------------------------------------------------#
# STEP 2: Output selection information
#------------------------------------------------------------------------------#
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
alp <- quantile(all_alpha, c(0.5, 0.025, 0.975), na.rm = T)
all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
h2 <- quantile(all_h2, c(0.5, 0.025, 0.975))
all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
p <- quantile(all_p, c(0.5, 0.025, 0.975))

df_t <- cbind(pheno, h2[1],h2[2],h2[3], p[1], p[2], p[3], alp[1], alp[2], alp[3])
colnames(df_t) <- c("pheno", "h2_est_1", "h2_est_2", "h2_est_3", "p_est_1", "p_est_2", "p_est_3", "alpha_est_1", "alpha_est_2", "alpha_est_3")
df_t <- as.data.frame(df_t)
write.table(df_t, out_alpha, sep='\t')
cat("Alphas written in :", out_alpha, "\n")


#------------------------------------------------------------------------------#
# STEP 3: Output estimated betas
#------------------------------------------------------------------------------#
sum_s_new = read.table(sum_s_p, hea=T)
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
stopifnot(length(beta_auto) == nrow(sum_s_new))

sum_s_new$effect_weight = beta_auto

gz_out = gzfile(out_beta, "wb")  # "wt" = write text
write.table(sum_s_new, gz_out, sep='\t', quote=F, row.names=F)
close(gz_out)
cat("Betas written in", out_beta, '\n')
