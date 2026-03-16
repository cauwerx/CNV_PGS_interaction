# Merges LoF per chromosome into single file and computes LoF autosomal burden

################################################################################
# Libraries
################################################################################
import pandas as pd

################################################################################
# STEP 1: Load LoF per CHR and concatenate
################################################################################

rare_burden = pd.read_csv('/mnt/project/penetrance_at_scale/data/24_countfiles_perCHR_MAF0.01/chr1.pLoF.sscore', sep='\t')
rare_burden.rename(columns={'NAMED_ALLELE_DOSAGE_SUM':'chr1'}, inplace=True)
rare_burden = rare_burden[['IID', 'chr1']]

for chr_num in range(2,23):
    chr_df = pd.read_csv('/mnt/project/penetrance_at_scale/data/24_countfiles_perCHR_MAF0.01/chr'+str(chr_num)+'.pLoF.sscore', sep='\t')
    chr_df.rename(columns={'NAMED_ALLELE_DOSAGE_SUM':'chr'+str(chr_num)}, inplace=True)
    rare_burden = pd.merge(rare_burden, chr_df[['IID', 'chr'+str(chr_num)]], on='IID', how='inner')

################################################################################
# STEP 2: Compute LoF burden across autosomes and save
################################################################################

rare_burden['LoF_total'] = rare_burden.drop(columns=['IID']).sum(axis=1)
rare_burden = rare_burden.sort_values(by='IID', ascending=True).reset_index(drop=True) # order by IID and re-index
rare_burden[['IID', 'LoF_total']].to_csv("LoF_burden.total.csv", header=True, index=False)
