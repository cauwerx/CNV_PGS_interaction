# Computes CNV effect size sign-concordance between 2 models: pheno~CNV and PGS~CNV
# Considers both pheno and PGS empirical correlation structure

#############################################################################################
# Libraries
#############################################################################################
install.packages("mvtnorm")
library(mvtnorm)
trials <- 1000
set.seed(84)

#############################################################################################
# STEP 1: Load empirical covariance matrices and effect sizes from each model
#############################################################################################
DIR <- "/mnt/project/sign_concordance/data/"

param_df <- read.csv(paste0(DIR, 'additive_model.parameters.csv'))
pheno_cov_matrix <- as.matrix(read.csv(paste0(DIR, 'pheno.cov_matrix.csv'), row.names=1))
pgs_cov_matrix <- as.matrix(read.csv(paste0(DIR, 'pgs.cov_matrix.csv'), row.names=1))

# Obtain PGS_mean from param_df, since we have same trait associated to different CNVs remove phenotype duplicates
# duplicated function identifies elements or rows that are duplicates of previous occurrences and extracts first occurence
PGS_mean_df <- param_df[!duplicated(param_df$phenotype), c('phenotype', 'PGS_mean')]

# Check if phenotypes are arranged in same order, unique function keeps order of appearance
stopifnot(identical(PGS_mean_df$phenotype, unique(param_df$phenotype)))
stopifnot(identical(PGS_mean_df$phenotype, rownames(pheno_cov_matrix)))
# Check matrices are same size and same number of unique phenotypes
stopifnot(nrow(pgs_cov_matrix) == nrow(pheno_cov_matrix))
stopifnot(nrow(pgs_cov_matrix) == nrow(PGS_mean_df))

pairs_n <- nrow(param_df) # 119 CNV-trait pairs
sample_n <- as.integer(mean(param_df$n))

#############################################################################################
# STEP 2: Simulate CNV, PGS and pheno. Then compute ratio of same-sign effect sizes.
#############################################################################################

# To store ratio of CNV->PGS and CNV->trait effects with same sign
random_ratios <- numeric(trials)

# Repeat for a given number of trials
for (i in 1:trials){
  if (i %% 100 == 0) print(i)
  c1_list <- numeric(pairs_n)
  c2_list <- numeric(pairs_n)

  # 1) Simulate PGS for all phenotypes
  # NOTE: PGS_mean_df$PGS_mean is a vector of PGS_mean for all pheno
  PGS_matrix  <- mvtnorm::rmvnorm(sample_n, mean=PGS_mean_df$PGS_mean, sigma=pgs_cov_matrix)
  colnames(PGS_matrix) <- PGS_mean_df$phenotype
  # 2) Simulate error term for all phenotypes
  e_matrix <- mvtnorm::rmvnorm(sample_n, mean=rep(0, nrow(pheno_cov_matrix)), sigma=pheno_cov_matrix)
  colnames(e_matrix) <- PGS_mean_df$phenotype

  # Iterate over 119 CNV-trait pairs
  for (j in 1:pairs_n){
    pheno_name <- param_df[j, 'phenotype']
    pheno_sd <- param_df[j, 'pheno_sd']
    carriers_n <- param_df[j, 'CNV_n']
    a <- param_df[j, 'CNV_coefficient']
    b <- param_df[j, 'PGS_coefficient']
   
    # 3) Simulate CNV (bernoulli)
    CNV <- rbinom(sample_n, 1, carriers_n/sample_n)
    
    # 4) Select PGS for that pheno
    PGS <- PGS_matrix[, pheno_name]
    
    if (!(pheno_sd^2 > (a^2*var(CNV) + b^2*var(PGS)))) {
      stop("error: pheno_var < explained_var")
    }
    
    # 5) Select error for that pheno and scale 
    e_sd <- sqrt(pheno_sd^2 - (a^2*var(CNV) + b^2*var(PGS)))
    e <- e_sd * e_matrix[, pheno_name]/sd(e_matrix[, pheno_name])

    # 6) Simulate phenotypes, with a and b taken from real values
    pheno <- a*CNV + b*PGS + e
    
    # 7) Get CNV -> PGS and CNV -> trait effects through linear regression
    model_1 <- lm(pheno~CNV)
    model_2 <- lm(PGS~CNV)
    c1_list[j] <- coef(model_1)[["CNV"]]
    c2_list[j] <- coef(model_2)[["CNV"]]
  }
  # 8) Assess the ratio showing sign agreement across all 119 pairs
  random_ratios[i] <- sum(c1_list * c2_list > 0)/pairs_n
}

write.csv(data.frame(random_ratios), 'samesign_simulated_ratios.csv', quote=FALSE, row.names=FALSE)
