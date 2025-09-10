# Get column names for the UKB sample QC file (ukb_sqc_v2.txt)

################################################################################
# Libraries
################################################################################

# Install libraries
devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE) # latest development version

# Load libraries
library(data.table)
library(ukbtools)

################################################################################
# STEP 1: Load data
################################################################################

sample_qc <- as.data.frame(fread("/mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb_sqc_v2.txt", header = F))

# Obtain names using ukb_gen_sqc_names from ukbtools
my_sqc_data <- ukb_gen_sqc_names(sample_qc) 

# Save names as dataframe
df <- as.data.frame(names(my_sqc_data))

################################################################################
# STEP 2: Save data
################################################################################

write.table(df, file = "ukb_sqc_v2.HEADER.txt", row.names = FALSE, col.names=FALSE, quote=FALSE)

