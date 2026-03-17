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

# 1. Genebass annotation file
genebass = pd.read_csv("/project/data/Genebass_variant_annotations.txt", sep="\t")
genebass.rename(columns={'gene': 'genename'}, inplace=True)
genebass.rename(columns={'annotation': 'annotation_genebass'}, inplace=True)
genebass.rename(columns={'AF': 'AF_genebass'}, inplace=True)

genebass['chr'] = genebass['chr'].str.replace('chr', '')
genebass['chr'] = genebass['chr'].str.replace('X', '23')
genebass['chr'] = genebass['chr'].str.replace('Y', '24')
genebass['variantID'] = genebass['chr']+":"+genebass['pos'].astype(str)+":"+genebass['REF']+":"+genebass['ALT']

# 2. SnpEff annotation file
snpeff = pd.read_csv("/project/data/Snpeff_variant_annotations.txt", sep=" ", names=["SNP", "gene_name", "annotation"])
snpeff.rename(columns={'gene_name': 'genename_ensemblID'}, inplace=True)
snpeff.rename(columns={'SNP': 'variantID'}, inplace=True)
snpeff.rename(columns={'annotation': 'annotation_snpeff'}, inplace=True)
pattern = r'(.+)\((ENSG\d+)\)'
# Extract the gene name and Ensembl ID
snpeff[['genename', 'ensemblID']] = snpeff['genename_ensemblID'].str.extract(pattern)

# 3. Merge both databases
combined = pd.merge(genebass[['variantID', 'genename', 'annotation_genebass', 'AF_genebass']], snpeff[['variantID', 'genename', 'annotation_snpeff']], on=['variantID','genename'], how='outer')

# 4. AAF info
dtype = {'ID': str, 'REF': str, 'MAJOR':str, 'ALT': str, 'MINOR':str, 'AAF': float, 'MAF': float, 'obs_allele_count':float, 'obs_allele_count_recoded':float}
maf_wb = pd.read_csv("/project/data/SNP_MAF_AAF.all_variants.nonrelated_whitebritish.txt", dtype=dtype)
maf_wb.rename(columns={'ID': 'variantID'}, inplace=True)
maf_wb.rename(columns={'AAF': 'AAF_wb'}, inplace=True)
maf_wb.rename(columns={'MAF': 'MAF_wb'}, inplace=True)
maf_wb.head()

# 5. Merge annotation of patogenicity with MAF annotation. Only keep variants with annotation (left merge)
anno_maf = pd.merge(combined, maf_wb[['variantID', 'REF', 'ALT', 'AAF_wb', 'MAF_wb', 'obs_allele_count']], on=["variantID"], how="left")
anno_maf = anno_maf[['variantID', 'REF', 'ALT', 'AF_genebass', 'AAF_wb', 'MAF_wb', 'genename', 'annotation_genebass', 'annotation_snpeff']]
anno_maf.to_csv('/project/data/Genebass_SnpEff_variant_annotations.with_MAF_unrelated_WB.txt')

################################################################################
# STEP 2: Filter variants
################################################################################

# Filter for variants within the AAF threshold (since pathogenicity refers to ALT allele) and pLoF annotated
sub_anno_maf = anno_maf.query('(AAF_wb < @MAF_threshold) & (annotation_genebass == "pLoF") & (annotation_snpeff == "LoF" | annotation_snpeff.isna())')

# get CHR and ALT allele
sub_anno_maf["CHR"] = sub_anno_maf["variantID"].str.split(":").str.get(0)
sub_anno_maf["ALT"] = sub_anno_maf["variantID"].str.split(":").str.get(3)

# Set weight=1 (same weight for all variants)
sub_anno_maf.loc[:,'weight'] = 1

chr_list = set(sub_anno_maf['CHR'])

# Check no duplicates of same variant
assert(sub_anno_maf[sub_anno_maf['variantID'].duplicated(keep=False)].sort_values('variantID').shape[0]==0)

# Create scorefile per chromosome
for chr_num in chr_list:
    chr_df = sub_anno_maf.query('CHR==@chr_num')
    chr_df[['variantID', 'ALT', 'weight']].to_csv('../data/24_scorefiles_perCHR_MAF0.01/chr'+str(chr_num)+'.'+pathogenicity+'.genebass_snpeff.correct.txt', sep=' ', header=False, index=False)

