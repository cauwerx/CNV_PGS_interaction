# Compute within-couple phenotypic correlation (in test set)

################################################################################
# Libraries
################################################################################
import pandas as pd
import numpy as np
from scipy import stats

project_DIR = "/mnt/project/"

################################################################################
# STEP 1: Load couples and phenotype data
################################################################################

# 1. Couples data (0=F, 1=M)
couples = pd.read_csv(project_DIR + '/data/hh_pairs_filter.jenny_sjaarda.txt', sep='\t')
couples = couples.astype('int32')

# Condition to swap: if member1 is male, swap members so the member1 is female
condition = couples['HOUSEHOLD_MEMBER1_sex'] == 1 # mask: is member1 male?
# If condition holds, then swap members. Both operations have to be done at same time (swapping M to member2 and F to member1)
# np.where chooses first option (swap), unless it's False for condition, in which case it chooses second option (don't swap)
couples['HOUSEHOLD_MEMBER1'], couples['HOUSEHOLD_MEMBER2'] = np.where(condition, couples['HOUSEHOLD_MEMBER2'], couples['HOUSEHOLD_MEMBER1']), np.where(condition, couples['HOUSEHOLD_MEMBER1'], couples['HOUSEHOLD_MEMBER2'])
# Now all first members are F, so correct sex column as well
couples['HOUSEHOLD_MEMBER1_sex'], couples['HOUSEHOLD_MEMBER2_sex'] = 0, 1

# 2. Phenotype categories
category_df = pd.read_csv(project_DIR + '/data/unique_phenotypes_annotated.43.txt', sep='\t')
category_order = ['Metabolic', 'Hematologic', 'Hepatic', 'Musculoskeletal', 'Cardiopulmonary',  'Renal', 'Neuropsychiatric', 'Endocrine']
category_df['CATEGORY'] = pd.Categorical(category_df['CATEGORY'], categories=category_order, ordered=True)
category_df = category_df.sort_values(by=["CATEGORY", "LABEL"], ascending=[True, True])
ordered_phenotypes = list(category_df['LABEL'])
ordered_phenotypes_abbrev = list(category_df['PHENO'])

# Create a dictionary mapping PHENO (keys) to LABEL (values)
rename_dict = pd.Series(category_df.LABEL.values, index=category_df.PHENO).to_dict()
inv_rename_dict = pd.Series(category_df.PHENO.values, index=category_df.LABEL).to_dict()

# 3. INT+cov-corrected phenotypes and PGSs
pheno = pd.read_csv(project_DIR+"/data/pheno/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.tsv.tar.gz", sep='\t') # test set
pgs = pd.read_csv(project_DIR+"/data/PGS.out_of_sample/PGS_continuous.tsv.gz", sep='\t')

pheno_F = pheno.copy()
pheno_F.rename(columns=rename_dict, inplace=True)
pheno_F = pheno_F[['IID']+ordered_phenotypes]
pheno_F.columns = [col + "_F" for col in pheno_F.columns]

pheno_M = pheno.copy()
pheno_M.rename(columns=rename_dict, inplace=True)
pheno_M = pheno_M[['IID']+ordered_phenotypes]
pheno_M.columns = [col + "_M" for col in pheno_M.columns]

print("Number of F in couples but not in test set:", len(set.difference(set(couples.HOUSEHOLD_MEMBER1), set(pheno_F.IID_F))))
print("Number of M in couples but not in test set:", len(set.difference(set(couples.HOUSEHOLD_MEMBER2), set(pheno_M.IID_M))))

couples_num = len(couples[couples['HOUSEHOLD_MEMBER1'].isin(pheno_F['IID_F']) & couples['HOUSEHOLD_MEMBER2'].isin(pheno_M['IID_M'])])
print("Number of couples in test set:", couples_num)

# 4. Add phenotype data to couples df
couples = pd.merge(couples, pheno_F, left_on='HOUSEHOLD_MEMBER1', right_on='IID_F')
couples = pd.merge(couples, pheno_M, left_on='HOUSEHOLD_MEMBER2', right_on='IID_M')
print("Number of couples in merged set:", len(couples))
# Drop unneeded columns
couples = couples.drop(['HOUSEHOLD_MEMBER1', 'HOUSEHOLD_MEMBER2'], axis=1)
# Rename sex cols
couples = couples.rename(columns={'HOUSEHOLD_MEMBER1_sex': 'sex_F', 'HOUSEHOLD_MEMBER2_sex': 'sex_M'})

# 5. Save data
couples[['IID_F', 'IID_M', 'sex_F', 'sex_M', 'HOUSE_ID', 'tar_group', 'kinship']].to_csv('couples.test.csv', index=False)

################################################################################
# STEP 2: phenotype_F - phenotype_M correlation
################################################################################

# 1. Correlation and p-value matrices
# Selecting column names to compute corr (omit IID column)
colsF = list(pheno_F.columns[1:])
colsM = list(pheno_M.columns[1:])
cols = colsF + colsM

# 2. Compute the correlation coefficients between all columns in the dataframe, then subset to those between F and M
# .corr already excludes NA/null values
corr_matrix = couples[cols].corr(method=lambda x, y: stats.pearsonr(x, y).statistic).loc[colsF,colsM]
# 3. Compute the p-values. Workaround with -identity matrix for the p-values is needed, because self-correlations are always set to p-value=1, so subtract identity matrix (see https://github.com/pandas-dev/pandas/issues/25726).
pvalue_matrix = couples[cols].corr(method=lambda x, y: stats.pearsonr(x, y).pvalue) - np.eye(len(couples[cols].columns))
pvalue_matrix = pvalue_matrix.loc[colsF,colsM]

# 4. Rename rows and columns (remove _F and _M, keeping in mind F=rows, M=columns)
corr_matrix.columns = corr_matrix.columns.str.replace('_M', '')
corr_matrix.index = corr_matrix.index.str.replace('_F', '')
pvalue_matrix.columns = pvalue_matrix.columns.str.replace('_M', '')
pvalue_matrix.index = pvalue_matrix.index.str.replace('_F', '')

# 5. Supplementary table
# Correlation
long_corr = corr_matrix.reset_index().melt(id_vars=corr_matrix.index.name or 'index')
long_corr.rename(columns={'index': 'PHENO_F', 'variable': 'PHENO_M', 'value': 'correlation'}, inplace=True)
long_corr.replace(inv_rename_dict, inplace=True)
# P-values
long_pval = pvalue_matrix.reset_index().melt(id_vars=pvalue_matrix.index.name or 'index')
long_pval.rename(columns={'index': 'PHENO_F', 'variable': 'PHENO_M', 'value': 'P'}, inplace=True)
long_pval.replace(inv_rename_dict, inplace=True)
# Save correlation + p-val file
long_df = pd.merge(long_corr, long_pval, on=['PHENO_F', 'PHENO_M'], how='inner')
long_df.to_csv('AM_matrix.tsv', sep='\t', index=False)

