# Find permutation of trait-CNV associations that maximizes trait-CNV p-values
# using Hungarian algorithm
setwd("/penetrance_at_scale/")

################################################################################
# Libraries 
################################################################################
library('clue')
library('data.table')

################################################################################
# STEP 1: Load data
################################################################################

# 1. CNV_GWAS triplets
# One row per original trait-PGS-CNV triplet.
# CNVs and traits are allowed to occur multiple times.

triplets <- read.table("data/CNV_GWAS_signals_ranking_v1.txt", sep='\t', header=TRUE)
triplets <- subset(triplets, SEX=='All' & TYPE=='continuous' & CHR!='X')
# There is one row with two top models, assign the most general
triplets$TOP_MODEL[triplets$TOP_MODEL == "M-DUP"] <- "M"
# CNV names
triplets$CNV <- paste(triplets$PHENO, 'CNV', triplets$CHR, triplets$CNVR_START, triplets$CNVR_STOP, sep='_')
triplets[c("CHR", "CNVR_START", "CNVR_STOP")] <- lapply(triplets[c("CHR", "CNVR_START", "CNVR_STOP")], FUN=as.integer)
triplets <- triplets[order(triplets$CHR, triplets$CNVR_START),]
# PGS names for PGS_GW
triplets$PGS <- paste(triplets$PHENO, 'PGS', sep='_')
# Only keep triplets: pheno-cnv-pgs
triplets <- triplets[, c('PHENO', 'CNV', 'PGS')]
n <- nrow(triplets)
stopifnot(n==119)

# 2. Matrix of pheno-cnv correlations
# rows = unique traits, cols = unique CNV encodings
p_m <- as.matrix(fread("data/07_sensitivity_analyses/01.trait_cnv.corr_pval.csv"), rownames="trait")
r_m <- as.matrix(fread("data/07_sensitivity_analyses/01.trait_cnv.corr_coef.csv"), rownames="trait")

# Check columns and rows are named the same
stopifnot(all(triplets$PHENO %in% rownames(p_m)))
stopifnot(all(triplets$CNV %in% colnames(p_m)))

################################################################################
# STEP 2: Build the 119 x 119 cost matrix
################################################################################

# rows = traits (can be repeated), cols = CNV encodings (can be repeated)
# To preserve multiplicity, even if the trait or CNV encoding are the same they appear as several rows, columns

cost_m <- matrix(NA_real_, nrow=n, ncol=n)

for (i in seq_len(n)){
  trait_i <- triplets$PHENO[i]
  for (j in seq_len(n)){
    CNV_j <- triplets$CNV[j]
    p <- p_m[trait_i, CNV_j] # Extract p-value for that trait-cnv combination
    # To circumvent the issue of having rows and cols with same name use indices 
    cost_m[i, j] <- p
  }
}

# We expect CNV_GWAS significant associations to have the lowest cost (p-value) 

# Some final checks
stopifnot(!anyNA(p_m))
stopifnot(!anyNA(cost_m))

################################################################################
# STEP 3: Use Hungarian algorithm to solve the assignment problem 
################################################################################

# Assign the CNVs to traits such that the total cost (the objective function = p-value) is maximal
# Since we want CNV-PGS pairings that are less likely to be associated
assignment <- solve_LSAP(cost_m, maximum = TRUE) 

# Convert the mapping to CNV indices 
assignment <- as.integer(assignment)

################################################################################
# STEP 4: Create null set
################################################################################

triplets$new_CNV <- triplets$CNV[assignment]
triplets$CNV_index <- seq_len(n)
triplets$new_CNV_index <- assignment

# Get the original and new cost for the triplets
triplets$cost <-  mapply(FUN=function(t, c) p_m[t, c], triplets$PHENO, triplets$CNV)
triplets$new_cost <-  mapply(FUN=function(t, c) p_m[t, c], triplets$PHENO, triplets$new_CNV)

# Save the original and new trait-CNV correlations
triplets$trait_CNV.r <-  mapply(FUN=function(t, c) r_m[t, c], triplets$PHENO, triplets$CNV)
triplets$trait_newCNV.r <-  mapply(FUN=function(t, c) r_m[t, c], triplets$PHENO, triplets$new_CNV)

################################################################################
# STEP 5: Compare total costs
################################################################################

objective_original <- sum(triplets$cost)
objective_null <- sum(triplets$new_cost)

cat("Original sum of p-vals:", objective_original, "\n")
cat("Maximum null sum of p-vals:", objective_null, "\n")
stopifnot(objective_null >= objective_original)

################################################################################
# STEP 6: Final checks
################################################################################

# Check multiplicities are preserved
stopifnot(identical(sort(triplets$CNV), sort(triplets$new_CNV)))
# Check no original assignments are preserved (this is not guaranteed though by Hungarian algo)
stopifnot(!any(triplets$CNV == triplets$new_CNV))

################################################################################
# STEP 7: Save new assignments
################################################################################

fwrite(triplets[,c('PHENO', 'PGS', 'CNV', 'new_CNV', 'trait_CNV.r', 'trait_newCNV.r')], "data/07_sensitivity_analyses/02.trait_cnv.empirical_and_null_sets.csv")
