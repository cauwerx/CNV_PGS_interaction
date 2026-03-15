# Assortative mating model (predicted vs observed PGS-CNV correlation)

################################################################################
# Libraries
################################################################################

packages <- c(
  "data.table",
  "dplyr",
  "glue",
  "vcmeta", # to get se.cor function which uses Bonett's formula: (1-r^2)/sqrt(N-3)
  "ggplot2",
  "ggrepel",
  "MetBrewer",
  "showtext",
  "R.utils"
)

# Install missing packages
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg)
  }
}

# Load libraries
lapply(packages, FUN=library, character.only = TRUE)

project_DIR <- "/mnt/project/"
mode <- "base" # base, couple_exclusion, sex_interact, location_coord, winners_curse

################################################################################
# Functions 
################################################################################

# Variance of product of independent random variables
# Exact formula: var(X1*X2*...*Xn) = prod(var(X_i) + E[X_i]^2) - prod(E[X_i]^2)

productSE <- function(correlations, SEs) {
  variances <- SEs^2
  
  # Compute (var(X_i) + (E[X_i])^2) for each variable then product
  first_terms <- variances + correlations^2    
  first_product <- prod(first_terms)
  
  # Compute E[X_i]^2 for each variable then product
  second_terms <- correlations^2
  second_product <- prod(second_terms)
  
  # Compute the SE of the product
  variance <- first_product - second_product
  SE <- sqrt(variance)
    
  return(SE)
}

# Confidence interval from SE
CI_from_SE <- function(correlation, SE, alpha=0.05){
  z_crit <- qnorm(1 - alpha/2)  # e.g., ~1.96 for 95% CI
  lo <- correlation - z_crit * SE
  hi <- correlation + z_crit * SE
  return(c(ci_lower = lo, ci_upper = hi))
}

################################################################################
# STEP 1: Load data
################################################################################

# 1. Phenotype categories
category_df <- read.table(
  file = file.path(project_DIR, "/data/unique_phenotypes_annotated.43.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)
# Sort by abbreviation
category_df <- category_df[order(category_df$PHENO), ]

# 2. CNV signals (CNVR associated to a complex trait)
CNVR <- read.table(
    file = file.path(project_DIR, "/data/CNV_GWAS_signals_ranking_v1.txt"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)
# Select only continuous phenotypes, not sex-specific + remove chromosome X
CNVR <- subset(CNVR, TYPE == "continuous" & SEX == "All" & CHR != 'X')
CNVR$CHR <- as.integer(CNVR$CHR)
CNVR$CNVR_START <- as.integer(CNVR$CNVR_START)
CNVR$CNVR_STOP <- as.integer(CNVR$CNVR_STOP)
CNVR$TOP_POS <- as.integer(CNVR$TOP_POS)
CNVR$REGION <- paste(CNVR$CHR, CNVR$CNVR_START, CNVR$CNVR_STOP, sep = "_")
CNVR$REGION_TRANS <- paste(CNVR$CHR, CNVR$CNVR_START-250000, CNVR$CNVR_STOP+250000, sep = "_")

# Order by phenotype
CNVR <- CNVR[order(CNVR$PHENO, CNVR$CHR, CNVR$CNVR_START), ]
# There is one top model which is both mirror (M) and duplication (D). Use M since it's more general.
both_model_i <- which(CNVR$TOP_MODEL == "M-DUP")
CNVR[both_model_i, "TOP_MODEL"] <- "M"
stopifnot(all(CNVR$TOP_MODEL != "M-DUP"))

stopifnot(nrow(CNVR) == 119)

if (mode=="winners_curse"){
    # new P-value threshold
    CNVR <- subset(CNVR, P<0.05/(11804*3))
    stopifnot(nrow(CNVR) == 98)
}

# Store phenotype names and their abbreviations
ordered_phenotypes_abbrev <- unique(CNVR$PHENO)
ordered_phenotypes <- category_df[category_df$PHENO %in% ordered_phenotypes_abbrev, 'LABEL']

# 3. Read genotype data (PGS & CNV) and phenotypes
if (mode=="location_coord"){
   pheno_pgs_cnv_path <- "/data/pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.pheno_location_correction.csv"
} else if (mode=="sex_interact"){
    pheno_pgs_cnv_path <- "/data/pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.pheno_sex_correction.csv"
} else {
   pheno_pgs_cnv_path <- "/data/pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.csv"
}
pheno_pgs_cnvburden_path <- "/data/pheno_pgs_cnv.autosomes.non_sex_specific.csv"

pheno_pgs_cnv <- as.data.frame(fread(file = file.path(project_DIR, pheno_pgs_cnv_path)))
pheno_pgs_cnvburden <- as.data.frame(fread(file = file.path(project_DIR, pheno_pgs_cnvburden_path)))

# 4. Read couples data
couples <- as.data.frame(fread(file = file.path(project_DIR, '/data/couples.test.csv')))

# 5. Get data of F and M separately and rename columns
pheno_F <- pheno_pgs_cnv %>% select(all_of(c("IID", ordered_phenotypes_abbrev)))
colnames(pheno_F) <- paste0(colnames(pheno_F), "_F")
pheno_M <- pheno_pgs_cnv %>% select(all_of(c("IID", ordered_phenotypes_abbrev)))
colnames(pheno_M) <- paste0(colnames(pheno_M), "_M")

# 6. Merge couples with pheno and genome data
new_couples <- merge(couples,     pheno_F, by = "IID_F",  all.x = TRUE)
new_couples <- merge(new_couples, pheno_M, by = "IID_M",  all.x = TRUE)
new_couples.IDs <- c(new_couples$IID_M, new_couples$IID_F)

print(paste("Number of couples in test set:", nrow(new_couples)))

if (mode=="couple_exclusion"){
    # 7. Remove couples from pheno_pgs_cnv
    pheno_pgs_cnv <- subset(pheno_pgs_cnv, !(IID %in% new_couples.IDs))
}

# 8. Load CNV burden
CNV_burden <- as.data.frame(
    fread(file = file.path(project_DIR, '/data/CNV_burden/CNV_burden.txt.gz'))
)

# 9. Merge pheno, PGS and CNV_burden data
pheno_pgs_cnvburden <- merge(pheno_pgs_cnvburden, CNV_burden, by = "IID",  all.x = TRUE)

# 10. Load cytogenetic bands
bands <- read.csv(paste0(project_DIR, '/data/helper_files/cytogenetic_bands.csv'), sep=";")
bands <- subset(bands, CHR!='X')
bands$CB_PHENO <- paste(bands$CB, '&', bands$PHENO)
bands$CHR <- as.integer(bands$CHR)
bands$POS <- as.integer(bands$POS)
bands <- bands[,c("PHENO", "CHR", "POS", "CB", "CB_PHENO")]
CNVR <- left_join(CNVR, bands, by = c("PHENO" = "PHENO", "CHR" = "CHR", "TOP_POS" = "POS"))

################################################################################
# STEP 2: Pheno-pheno couples correlation
################################################################################

# 1. Create an empty data frame for corr, CI, p-values & n
corr_df.v1 <- data.frame(matrix(NA, nrow = length(ordered_phenotypes_abbrev), ncol = 0), row.names = ordered_phenotypes_abbrev)

# 2. Add pheno_abbrev as another column, for labeling purposes
corr_df.v1$pheno_abbrev <- rownames(corr_df.v1)
corr_df.v1$pheno_abbrev <- gsub("_", " ", corr_df.v1$pheno_abbrev)
# Capitalize first letter, if it's not already capital
corr_df.v1$pheno_abbrev <- ifelse(grepl("^[a-z]", corr_df.v1$pheno_abbrev),
                 sub("^(.)", "\\U\\1", corr_df.v1$pheno_abbrev, perl = TRUE),
                 corr_df.v1$pheno_abbrev)

# 3. Add phenotype names as another column, for labeling purposes
corr_df.v1$pheno_name <- ordered_phenotypes

# 4. Iterate through each phenotype and compute correlation + p-values
for (pheno in ordered_phenotypes_abbrev) {
    # For 'pheno-pheno'
    # Select the corresponding phenotype columns and only non-NaN values
    valid_data <- new_couples %>%
        select(all_of(c(paste0(pheno, "_F"), paste0(pheno, "_M")))) %>%
        na.omit()

    if (nrow(valid_data) > 0) {
        # Compute correlations for 'pheno_F-pheno_M'. 1=F, 2=M
        test <- cor.test(valid_data[[1]], valid_data[[2]], method = "pearson")
        corr_df.v1[pheno, "pheno_F-pheno_M"] <- test$estimate
        corr_df.v1[pheno, "pheno_F-pheno_M.SE"] <- se.cor(test$estimate, s=0, n=nrow(valid_data))[,2]
        corr_df.v1[pheno, "pheno_F-pheno_M.CI1"] <- test$conf.int[[1]]
        corr_df.v1[pheno, "pheno_F-pheno_M.CI2"] <- test$conf.int[[2]]
        corr_df.v1[pheno, "pheno_F-pheno_M.p"] <- test$p.value
        corr_df.v1[pheno, "pheno_F-pheno_M.n"] <- nrow(valid_data)
    }
}

################################################################################
# STEP 3: Total vs Indirect CNV-PGS correlation (119 CNV-trait pairs)
################################################################################

# 1. Computing total and indirect pgs-cnv correlations
# create dataframe to save results for CNV-trait pairs
cnv_pgs.corr.v1 <- data.frame()

for (i in seq_len(nrow(CNVR))) {
  row_i <- CNVR[i,]
  pheno <- row_i[["PHENO"]]
  region <- row_i[["REGION"]]
  region_trans <- row_i[["REGION_TRANS"]]
  CB_PHENO <- row_i[["CB_PHENO"]]
  CB <- row_i[["CB"]]

  # Get pheno-PGS correlation info
  pgs_name <- paste0(pheno, '_chr', region_trans, '_PGS')
  valid_data <- subset(pheno_pgs_cnv, !is.na(pheno_pgs_cnv[[pheno]]) & !is.na(pheno_pgs_cnv[[pgs_name]]))
  if (nrow(valid_data) > 0) {
    test <- cor.test(valid_data[[pheno]], valid_data[[pgs_name]], method = "pearson")
    r1 <- unname(test$estimate)
    p1 <- test$p.value
    n1 <- nrow(valid_data)
    se1 <- se.cor(r1, s=0, n=n1)[,2]
  }

  # Get pheno_F-pheno_M correlation info
  r2  <- corr_df.v1[pheno, "pheno_F-pheno_M"]
  p2  <- corr_df.v1[pheno, "pheno_F-pheno_M.p"]
  n2  <- corr_df.v1[pheno, "pheno_F-pheno_M.n"]
  se2 <- corr_df.v1[pheno, "pheno_F-pheno_M.SE"]

  # Compute pheno-cnv correlation
  cnv_status <- paste0(pheno, '_CNV_', region)
  valid_data <- subset(pheno_pgs_cnv, !is.na(pheno_pgs_cnv[[pheno]]) & !is.na(pheno_pgs_cnv[[cnv_status]]))
  if (nrow(valid_data) > 0) {
    test <- cor.test(valid_data[[pheno]], valid_data[[cnv_status]], method = "pearson")
    r3 <- unname(test$estimate)
    p3 <- test$p.value
    n3 <- nrow(valid_data)
    se3 <- se.cor(r3, s=0, n=n3)[,2]
  }

  # Compute pgs-cnv correlation (total)
  valid_data <- subset(pheno_pgs_cnv, !is.na(pheno_pgs_cnv[[pgs_name]]) & !is.na(pheno_pgs_cnv[[cnv_status]]))
  if (nrow(valid_data) > 0) {
    test <- cor.test(valid_data[[pgs_name]], valid_data[[cnv_status]], method = "pearson")
    r4 <- unname(test$estimate)
    p4 <- test$p.value
    n4 <- nrow(valid_data)
    se4 <- se.cor(r4, s=0, n=n4)[,2]
    r4_lo <- test$conf.int[[1]]
    r4_up <- test$conf.int[[2]]
  }

   # Compute pgs-cnv correlation (indirect)
    r_indirect <- r1*r2*r3
    indirect_SE <- productSE(c(r1, r2, r3), c(se1, se2, se3))
    indirect_CI <- CI_from_SE(r_indirect, indirect_SE)

    # pheno without _
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno_abbrev"] <- gsub("_", " ", pheno)
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "PHENO"] <- pheno
    # bands info
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "CB_PHENO"] <- CB_PHENO
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "CB"] <- CB
    # gamma
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-PGS"] <- r1
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-PGS.SE"] <- se1
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-PGS.p"] <- p1
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-PGS.n"] <- n1
    # beta
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno_F-pheno_M"] <- r2
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno_F-pheno_M.SE"] <- se2
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno_F-pheno_M.p"] <- p2
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno_F-pheno_M.n"] <- n2
    # alpha
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-CNV"] <- r3
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-CNV.SE"] <- se3
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-CNV.p"] <- p3
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "pheno-CNV.n"] <- n3
    # total r(CNV, PGS)
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "total"] <- r4
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "total.SE"] <- se4
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "total.p"] <- p4
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "total.n"] <- n4
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "total.CI1"] <- r4_lo
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "total.CI2"] <- r4_up
    # indirect r(CNV, PGS)
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "indirect"] <- r_indirect
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "indirect.SE"] <- indirect_SE
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "indirect.CI1"] <- indirect_CI[['ci_lower']]
    cnv_pgs.corr.v1[paste0(pheno, "_CNV_", region), "indirect.CI2"] <- indirect_CI[['ci_upper']]

}

fwrite(cnv_pgs.corr.v1[, !(names(cnv_pgs.corr.v1) %in% c("pheno_abbrev", "CB_PHENO"))],
       paste0("path_coefficients.cnv_trait_pairs.pgs_trans.", mode, ".tsv"), col.names = T, row.names = T, quote = F, sep = "\t")


if (mode=="base"){

################################################################################
# STEP 4: Total vs Indirect CNV burden-PGS correlation
################################################################################

# 1. Load phenotypes significantly affected by burden
pheno_affected <- c("neuroticism", "GS", "vitamin_D", "WHR", "fluid_intelligence", "WHRadjBMI", "heart_rate", "neutrophil_count", "FVC",  "body_fat_mass", "WBC_count", "BMI", "CRP", "ApoA", "HbA1c", "cystatinC", "BMD", "HDL", "ALP", "platelet_count", "height")

# 2. Load pheno ~ CNV_burden model
CNVburden_model <- as.data.frame(fread(file.path(project_DIR, "/55_out_of_sample.analyses/data/CNVBurden_model.txt"), header=TRUE))
# Check pheno_affected are the significant ones:
stopifnot(length(setdiff(subset(CNVburden_model,P_CNV_GENES<0.05/43)$PHENO, pheno_affected))==0)

# 1. Computing total and indirect pgs-cnv_burden correlations
# Create dataframe to save results
cnv_pgs.corr.v2 <- data.frame()
cnv_burden <- "BURDEN_GENES"

for (pheno in pheno_affected) {  
  pgs_name <- paste0(pheno, "_PGS")
    
  # Get pheno-PGS correlation info
  valid_data <- subset(pheno_pgs_cnvburden, !is.na(pheno_pgs_cnvburden[[pheno]]) & !is.na(pheno_pgs_cnvburden[[pgs_name]]))
  if (nrow(valid_data) > 0) {
    test <- cor.test(valid_data[[pheno]], valid_data[[pgs_name]], method = "pearson")  
    r1 <- unname(test$estimate)
    p1 <- test$p.value
    n1 <- nrow(valid_data)
    se1 <- se.cor(r1, s=0, n=n1)[,2]
  }

  # Get pheno_F-pheno_M correlation info
  valid_data <- new_couples %>%
        select(all_of(c(paste0(pheno, "_F"), paste0(pheno, "_M")))) %>%
        na.omit()
    
  if (nrow(valid_data) > 0) {
        # Compute correlations for 'pheno_F-pheno_M'. 1=F, 2=M
        test <- cor.test(valid_data[[1]], valid_data[[2]], method = "pearson")
        r2 <- unname(test$estimate)
        p2 <- test$p.value
        n2 <- nrow(valid_data)
        se2 <- se.cor(r2, s=0, n=n2)[,2]
    }
    
  # Compute pheno-cnv_burden correlation
  valid_data <- subset(pheno_pgs_cnvburden, !is.na(pheno_pgs_cnvburden[[pheno]]) & !is.na(pheno_pgs_cnvburden[[cnv_burden]]))
  if (nrow(valid_data) > 0) {
    test <- cor.test(valid_data[[pheno]], valid_data[[cnv_burden]], method = "pearson")  
    r3 <- unname(test$estimate)
    p3 <- test$p.value
    n3 <- nrow(valid_data)
    se3 <- se.cor(r3, s=0, n=n3)[,2]
  }
      
  # Compute pgs-cnv_burden correlation (total) 
  valid_data <- subset(pheno_pgs_cnvburden, !is.na(pheno_pgs_cnvburden[[pgs_name]]) & !is.na(pheno_pgs_cnvburden[[cnv_burden]]))
  if (nrow(valid_data) > 0) {
    test <- cor.test(valid_data[[pgs_name]], valid_data[[cnv_burden]], method = "pearson")
    r4 <- unname(test$estimate)
    p4 <- test$p.value
    n4 <- nrow(valid_data)
    se4 <- se.cor(r4, s=0, n=n4)[,2]
    r4_lo <- test$conf.int[[1]]
    r4_up <- test$conf.int[[2]]      
  }    
    
    # Compute pgs-cnv_burden correlation (indirect) 
    r_indirect <- r1*r2*r3    
    indirect_SE <- productSE(c(r1, r2, r3),c(se1, se2, se3))
    indirect_CI <- CI_from_SE(r_indirect, indirect_SE)
    
    # pheno without _        
    cnv_pgs.corr.v2[pheno, "pheno_abbrev"] <- gsub("_", " ", pheno)  
    cnv_pgs.corr.v2[pheno, "PHENO"] <- pheno

    # gamma 
    cnv_pgs.corr.v2[pheno, "pheno-PGS"] <- r1
    cnv_pgs.corr.v2[pheno, "pheno-PGS.SE"] <- se1
    cnv_pgs.corr.v2[pheno, "pheno-PGS.p"] <- p1
    cnv_pgs.corr.v2[pheno, "pheno-PGS.n"] <- n1
    # beta 
    cnv_pgs.corr.v2[pheno, "pheno_F-pheno_M"] <- r2
    cnv_pgs.corr.v2[pheno, "pheno_F-pheno_M.SE"] <- se2
    cnv_pgs.corr.v2[pheno, "pheno_F-pheno_M.p"] <- p2
    cnv_pgs.corr.v2[pheno, "pheno_F-pheno_M.n"] <- n2
    # alpha     
    cnv_pgs.corr.v2[pheno, "pheno-CNVburden"] <- r3
    cnv_pgs.corr.v2[pheno, "pheno-CNVburden.SE"] <- se3
    cnv_pgs.corr.v2[pheno, "pheno-CNVburden.p"] <- p3
    cnv_pgs.corr.v2[pheno, "pheno-CNVburden.n"] <- n3
    # total r(CNV, PGS) 
    cnv_pgs.corr.v2[pheno, "total"] <- r4
    cnv_pgs.corr.v2[pheno, "total.SE"] <- se4
    cnv_pgs.corr.v2[pheno, "total.p"] <- p4
    cnv_pgs.corr.v2[pheno, "total.n"] <- n4
    cnv_pgs.corr.v2[pheno, "total.CI1"] <- r4_lo
    cnv_pgs.corr.v2[pheno, "total.CI2"] <- r4_up    
    # indirect r(CNV, PGS)     
    cnv_pgs.corr.v2[pheno, "indirect"] <- r_indirect
    cnv_pgs.corr.v2[pheno, "indirect.SE"] <- indirect_SE
    cnv_pgs.corr.v2[pheno, "indirect.CI1"] <- indirect_CI[['ci_lower']]
    cnv_pgs.corr.v2[pheno, "indirect.CI2"] <- indirect_CI[['ci_upper']]
    
}
   
fwrite(cnv_pgs.corr.v2[, !(names(cnv_pgs.corr.v2) %in% c("pheno_abbrev"))], "path_coefficients.cnv_burden.tsv", col.names = T, row.names = T, quote = F, sep = "\t") 

}

