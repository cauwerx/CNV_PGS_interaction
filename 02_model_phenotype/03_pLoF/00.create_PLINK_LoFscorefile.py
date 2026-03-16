# Create scorefiles for PLINK

################################################################################
# Libraries
################################################################################
import pandas as pd

################################################################################
# Parameters
################################################################################
MAF_threshold = 0.01
pathogenicity = "pLoF" 

################################################################################
# STEP 1: Load data
################################################################################

# 1. Genebass + SnpEff annotation file
anno_path = "./data/SnpEff_Genebass."+pathogenicity+".with_MAF.nonrelated_whitebritish.txt"
dtype = {'SNP': str, 'genename': str, 'annotation_snpeff': str, 'annotation_genebass': str, 'MAF_wb': float, 'AAF_wb': float}
anno = pd.read_csv(anno_path, dtype=dtype)

# get CHR, ALT, Gene and weight=1 (same weight for all variants)
anno["CHR"] = anno["SNP"].str.split(":").str.get(0).astype(int)
anno["ALT"] = anno["SNP"].str.split(":").str.get(3)

################################################################################
# STEP 2: Filter variants
################################################################################

# Filter for variants within the MAF threshold annotated
sub_anno = anno.query('(MAF_wb < @MAF_threshold)')

# Check no duplicate variants
assert(sub_anno[sub_anno['SNP'].duplicated(keep=False)].shape[0]==0)

# Set weight=1 (same weight for all variants)
sub_anno.loc[:,'weight'] = 1

chr_list = set(sub_anno['CHR'])

# Create scorefile per chromosome
for chr_num in chr_list:
    chr_df = sub_anno.query('CHR==@chr_num')
    chr_df[['SNP', 'ALT', 'weight']].to_csv('./data/chr'+str(chr_num)+'.'+pathogenicity+'.genebass_snpeff.txt', sep=' ', header=False, index=False)

