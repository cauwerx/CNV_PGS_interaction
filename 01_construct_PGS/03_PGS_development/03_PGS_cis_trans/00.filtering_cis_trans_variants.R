################################################################################
# Libraries
################################################################################
install.packages('readxl')
library('readxl')

################################################################################
# Parameters
################################################################################
setwd('/project/data/')
type <- "250kb" # 0kb, 50kb, 100kb, 250kb, CHR
out_DIR <- './PGS_scorefiles.cis_trans/'
output_file <- paste0(out_DIR, "/_out_files/PGS_partitioning_cnv_signals_", type ,".out") 
# Create output file and then create connection to it for appending
writeLines(paste("Type of filtering:", type), output_file)
fileConn<-file(output_file, "a")

################################################################################
# STEP 1: Define cis/trans CNVR for each CNV-trait pair 
################################################################################

# 1.1) Load CNV-GWAS signals
CNVR <- read.table("./CNV_GWAS_signals_ranking_v1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 1.2) Filter for continuous phenotypes, only autosomal CHR
CNVR <- subset(CNVR, TYPE == "continuous" & CHR != 'X')
CNVR$CHR <- as.integer(CNVR$CHR)
CNVR$CNVR_START <- as.integer(CNVR$CNVR_START)
CNVR$CNVR_STOP <- as.integer(CNVR$CNVR_STOP)

# 1.3) Obtain CNVR limits (REGION_ID)
if (type %in% c("0kb","50kb", "100kb", "250kb")) {
  # Get size
  size_kb <- as.numeric(gsub("kb", "", type))
  # Add flanking window
  CNVR$CNVR_START <- CNVR$CNVR_START-size_kb*1000
  CNVR$CNVR_STOP <- CNVR$CNVR_STOP+size_kb*1000
  CNVR$REGION_ID <- paste(CNVR$PHENO, paste0('chr', CNVR$CHR), CNVR$CNVR_START, CNVR$CNVR_STOP, sep = "_")
  CNVR <- CNVR[order(CNVR$PHENO, CNVR$CHR, CNVR$CNVR_START), ]
  cnv_signals_df <- CNVR[,c('PHENO', 'CHR', 'CNVR_START', 'CNVR_STOP', 'REGION_ID')]
} else if (type == "CHR") {
  CNVR$REGION_ID <- paste(CNVR$PHENO, 'chr', sep = "_")
  CNVR <- CNVR[order(CNVR$PHENO, CNVR$CHR), ]
  cnv_signals_df <- unique(CNVR[,c('PHENO', 'CHR', 'REGION_ID')])
}

n_pheno <- length(unique(cnv_signals_df$PHENO))
writeLines(paste(n_pheno, 'unique continuous phenotypes with CNV signal'), fileConn)
writeLines(paste(nrow(cnv_signals_df), 'unique signals\n'), fileConn)

# 1.4) Save this file with REGION_ID
write.table(cnv_signals_df, paste0(out_DIR, "/_out_files/PGS_partitioning_cnv_signals_", type ,".txt"), sep="\t", quote = FALSE, row.names=FALSE)

################################################################################
# STEP 2: Load PGS files and separate headers from scores
################################################################################

# 2.1) Read the corresponding PGS files, one PGS per phenotype
PGS_list <- unique(cnv_signals_df$PHENO)

# 2.2) Initialize list of dataframes 
#      - names = phenotype
#      - values = header of scorefiles (PGS_header_dict)
#               = weights for all SNPs (PGS_dict)
#               = weights for cis SNPs (PGS_dict_in) 
#               = weights for trans SNPs (PGS_dict_out)

# input
PGS_header_dict <- list()
PGS_dict <- list()
# output
PGS_dict_in <- list()
PGS_dict_out <- list()

# 2.3) Populate PGS_header_dict and PGS_dict
for (pheno in PGS_list) {
  file_path <- paste0("./PGS_scorefiles.parsed/", pheno, ".betas.tsv")
  # Read the header and the table in each scorefile
  PGS_header_dict[[pheno]] <-  readLines(file_path, n = 4)
  PGS_dict[[pheno]] <- read.table(file_path, skip=4, sep='\t', header=TRUE)
  PGS_header_dict[[pheno]][5] <- paste0("#variants_number=", nrow(PGS_dict[[pheno]]))
}

################################################################################
# STEP 3: Define PGS cis/trans scorefiles
################################################################################

# 3.1) Iterate over phenotype, then over CNVR
for (pheno in PGS_list){
  print(pheno)
  writeLines("-----------------------------------------------------------------------------------", fileConn)
  writeLines("", fileConn)
  writeLines(paste("Processing phenotype:", pheno), fileConn)
  # Number of variants in PGS
  var_n <- nrow(PGS_dict[[pheno]])
  writeLines(paste("Number of variants in PGS:", var_n), fileConn)
  writeLines("", fileConn)
  # Number of lines in header
  header_n <- length(PGS_header_dict[[pheno]])
  # Get the chr numbers where the phenotype-associated CNV appears (a CNV can be in several regions and associated to same phenotype)
  cnv_regions <- subset(cnv_signals_df, PHENO==pheno)
  writeLines(paste0("Regions associated to ", pheno, ":"), fileConn)
  write.table(cnv_regions, fileConn, sep="\t", quote = FALSE, append=TRUE, row.names=FALSE)
  writeLines("", fileConn)
  chromosomes <- unique(cnv_regions$CHR)
  # Iterate over each unique CHR where CNV is located 
  for (num in chromosomes){
    # Select rows for that particular CHR
    df <- subset(cnv_regions, CHR==num)
    writeLines(paste0("Regions associated to ", pheno, " in chr ", num, ":"), fileConn)
    write.table(df, fileConn, sep="\t", quote = FALSE, append=TRUE, row.names=FALSE)
    writeLines("", fileConn)
    if (type %in% c("0kb","50kb", "100kb", "250kb")){
      # For all CNVs in that particular CHR (iterate over all CNVs in a single CHR)
      for (i in (1:nrow(df))){
        start <- df[i,]$CNVR_START
        end <- df[i,]$CNVR_STOP
        cnv_positions <- c(start:end)
        writeLines(paste("Range:", start, "-",end), fileConn)
        # in SNPs = those in that chr and that range
        # out SNPS = those outside that chr, or inside that chr but outside that range
        in_df <- subset(PGS_dict[[pheno]], chr_name==num & chr_position %in% cnv_positions)
        out_df <- subset(PGS_dict[[pheno]], (chr_name==num & !(chr_position %in% cnv_positions)) | chr_name!=num)
        in_file_name <- paste0(out_DIR, df$REGION_ID[i], "_", type, "_", "in.txt")
        out_file_name <- paste0(out_DIR, df$REGION_ID[i], "_", type, "_", "out.txt")
        # Change number of variants in the header of in_file and save header + in_df
        in_var_n <- nrow(in_df)
        PGS_header_dict[[pheno]][header_n] <- paste0("#variants_number=", in_var_n)
        writeLines(PGS_header_dict[[pheno]],  in_file_name)
        write.table(in_df, in_file_name, sep="\t", quote = FALSE, append=TRUE, row.names=FALSE)
        writeLines(paste("Saving:", in_file_name), fileConn)
        # Change number of variants in the header of out_file and save header + out_df
        out_var_n <- nrow(out_df)
        PGS_header_dict[[pheno]][header_n] <- paste0("#variants_number=", out_var_n)
        writeLines(PGS_header_dict[[pheno]],  out_file_name)
        write.table(out_df, out_file_name, sep="\t", quote = FALSE, append=TRUE, row.names=FALSE)
        writeLines(paste("Saving:", out_file_name), fileConn)
        # Check variants add up to total
        stopifnot(in_var_n + out_var_n == var_n)
        writeLines(paste0("IN: ", in_var_n, ", OUT: ", out_var_n, ", TOTAL: ", var_n), fileConn)
        writeLines("", fileConn)
      }
    } else if (type=="CHR") {
      # in SNPs = those in that chr
      # out SNPS = those outside that chr
      in_df <- subset(PGS_dict[[pheno]], chr_name==num)
      out_df <- subset(PGS_dict[[pheno]], chr_name!=num)
      in_file_name <- paste0(out_DIR, df$REGION_ID, df$CHR, "_", "in.txt")
      out_file_name <- paste0(out_DIR, df$REGION_ID, df$CHR, "_", "out.txt")
      # Change number of variants in the header of in_file and save header + in_df
      in_var_n <- nrow(in_df)
      PGS_header_dict[[pheno]][header_n] <- paste0("#variants_number=", in_var_n)
      writeLines(PGS_header_dict[[pheno]],  in_file_name)
      write.table(in_df, in_file_name, sep="\t", quote = FALSE, append=TRUE, row.names=FALSE)
      writeLines(paste("Saving:", in_file_name), fileConn)
      # Change number of variants in the header of out_file and save header + out_df
      out_var_n <- nrow(out_df)
      PGS_header_dict[[pheno]][header_n] <- paste0("#variants_number=", out_var_n)
      writeLines(PGS_header_dict[[pheno]],  out_file_name)
      write.table(out_df, out_file_name, sep="\t", quote = FALSE, append=TRUE, row.names=FALSE)
      writeLines(paste("Saving:", out_file_name), fileConn)
      # Check variants add up to total
      stopifnot(in_var_n + out_var_n == var_n)
      writeLines(paste0("IN: ", in_var_n, ", OUT: ", out_var_n, ", TOTAL: ", var_n), fileConn)
      writeLines("", fileConn)
    }
  } 
}

close(fileConn)
