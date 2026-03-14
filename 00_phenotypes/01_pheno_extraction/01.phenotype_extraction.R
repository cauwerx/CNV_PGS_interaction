############################################################################################################
# Libraries
############################################################################################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)

############################################################################################################
# STEP 1: Load Main Data 
############################################################################################################
DIR="./project/"
OUT_DIR="./project/data/04_phenotypes/"

# File with columns: FieldID, phenotype_abbreviation
pheno_list <- as.data.frame(fread(paste0(DIR,"/data/continuous_traits.csv"), header = T))


############################################################################################################
# STEP 2: Order raw phenotype files 
############################################################################################################

# full.names includes path, recursive includes subdirectories. it returns a df of files and metadata
raw_pheno_file <- file.info(list.files(paste0(DIR,"/uk_biobank/phenotypes/"), pattern = "^ukb[0-9]+.*\\.csv", full.names = T, recursive = F))
# Save only path and date
raw_pheno_file <- data.frame(File = rownames(raw_pheno_file), Date = raw_pheno_file[,4])
raw_pheno_file <- raw_pheno_file[order(raw_pheno_file$Date, decreasing = T), ]
print(paste0("There are ", nrow(raw_pheno_file), " raw phenotype files"))

############################################################################################################
# STEP 3: Identify most recent location 
############################################################################################################

pheno_list$File <- NA
pheno_list$Col <- NA

# Iterate over each phenotype
for (p in 1:nrow(pheno_list)) {
        # Define pheno
        pheno <- pheno_list[p, "Pheno"]
        ID <- pheno_list[p, "FieldID"]
        print(paste0("Extracting most recent file for ", pheno))
	    # Counter to iterate over files in raw_pheno_file. Start always from first file for each pheno.
        counter <- 1

        # Define most recent location --> loop through raw files
	    # File is NA until it's found in a file (first match ends the while loop, but this file is the most recent one since they're ordered by creation time) 
        while (is.na(pheno_list[p, "File"])) {
		ukb_filename <- as.character(raw_pheno_file[counter, "File"])
		print(paste0('Looking at file:', ukb_filename))
		        # Don't read any rows, but still parse column names
                header <- fread(ukb_filename, nrow = 0)
		        # Search for phenotype "<ID>-" among the list of column_names, and return indices in such list
                col <- grep(paste0("^", ID, "-"), names(header))
                if (length(col) > 0) {
			            # Save name of file where it was found
                        pheno_list[p, "File"] <- ukb_filename
			            # Print only basename
                        print(paste0("Most recent location: ", sub(".*/", "", pheno_list[p, "File"])))
			            # Save indices (# of the columns where it was found)
                        pheno_list[p, "Col"] <- paste(col, collapse = "_")}
                counter <- counter +1
        }
}
rm(p, pheno, ID, counter, header, col)
fwrite(pheno_list, paste0(OUT_DIR, "pheno_list.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

############################################################################################################
# STEP 4: Extract phenotypes + Average 
############################################################################################################

# Read the eids only (second column), drop=F is to avoid getting a vector instead of df
phenotypes <- as.data.frame(fread(paste0(DIR,"/cauwerx/cauwerx/projects/general_data/UKBB/eid_no_redacted3.txt.gz"), header = T, select = c(2)), drop = F)

for(p in 1:nrow(pheno_list)) {

        # Define pheno
        pheno <- pheno_list[p, "Pheno"]
        ID <- pheno_list[p, "FieldID"]
        file <- pheno_list[p, "File"]
        col <- pheno_list[p, "Col"]
        print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))

        # Extract pheno columns
	    # Select the eid (first column) and all the phenotype measurements available
        temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))

        # Average (if more than one column present)
        if (ncol(temp_pheno) > 2) {
        temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))}

        # Merge current pheno with previous
        colnames(temp_pheno) <- c("eid", pheno)
        phenotypes <- full_join(phenotypes, temp_pheno, by = "eid")

}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))


############################################################################################################
# STEP 5: Composite traits
############################################################################################################

# Grip strength: mean left and right
phenotypes$GS <- (phenotypes$grip_strength_left + phenotypes$grip_strength_right)/2
phenotypes <- phenotypes[, !names(phenotypes) %in% c("grip_strength_left", "grip_strength_right")]

# WHR: waist circumference divided by hip circumference
phenotypes$WHR <- phenotypes$waist_circumference/phenotypes$hip_circumference
phenotypes <- phenotypes[, !names(phenotypes) %in% c("waist_circumference", "hip_circumference")]

# WHRadjBMI: waist circumference divided by hip circumference adjusted for BMI --> Done after selecting samples

print(paste0("Dimensions of phenotype table after adding composite traits: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))


############################################################################################################
# STEP 6: Replace non-conclusive entries by NA 
############################################################################################################

# -1: do not know; -2: only had twins; -3: prefer not to answer
phenotypes[phenotypes < 0] <- NA


############################################################################################################
# STEP 7: Select samples & add WHRadjBMI
############################################################################################################

# Sample QC
sample_qc <- as.data.frame(fread(paste0(DIR,"/uk_biobank/phenotypes/ukb_sqc_v2.txt" ), header = F, select = c(3:68)))
colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")

# Sample eid
sample_eid <- as.data.frame(fread(paste0(DIR, "/uk_biobank/genotypes/plink/_001_ukb_cal_chr1_v2.fam"), header = F, select = c(1,5), col.names = c("eid", "sex")))

# Merge eid and sampleQC Data
df <- cbind(sample_eid, sample_qc)
print(paste0("Start: ", nrow(df), " individuals"))

# Exclude non-white, non-British ancestry samples
df <- df[which(df$white_british == 1), ]
print(paste0("Exclude non-white, non-British ancestry samples: ", nrow(df), " individuals"))

# Exclude samples with non-matching submitted vs. inferred sex
df <- df[which(df$submitted_gender == df$inferred_gender), ]
print(paste0("Exclude sex mismatches: ", nrow(df), " individuals"))

pheno_wb <- merge(phenotypes, df[, c("eid", "sex")], by = "eid")
pheno_wb$WHRadjBMI <- residuals(lm(pheno_wb$WHR ~ pheno_wb$BMI, na.action = na.exclude))

colnames(pheno_wb)[1] <- "IID"
# re-arrange so that sex is before phenotypes
cols <- names(pheno_wb)
new_order <- c("IID", "sex", setdiff(cols, c("IID", "sex")))
pheno_wb <- pheno_wb[, new_order]

############################################################################################################
# STEP 8: Save file for all WB individuals
############################################################################################################

print(paste0("Dimensions of phenotype table for all white British individuals: ", ncol(pheno_wb), " pheno x ", nrow(pheno_wb), " eids"))
fwrite(pheno_wb, paste0(OUT_DIR, "pheno_continuous_WB_raw.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

