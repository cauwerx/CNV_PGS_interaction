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
pheno_list <- as.data.frame(fread("./coord_list.csv", header = T))

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
		print(pheno_list[p,])
        }
}
rm(p, pheno, ID, counter, header, col)
fwrite(pheno_list, paste0(OUT_DIR, "location_coord.list.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

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

fwrite(phenotypes, paste0(OUT_DIR, "location_coord.tsv"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


