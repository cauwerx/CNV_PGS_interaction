# Compare within-couple phenotypic correlation (original couple definition vs stricter couple definition)

##################################################################################################################
# Libraries
##################################################################################################################
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import norm
  
project_DIR = "/mnt/project/"

##################################################################################################################
# Functions 
##################################################################################################################

def pearson_se(r, n):
    """Standard error of Pearson r using vcmeta::se.cor R function, 
    corresponding to Bonett's formula: (1-r^2)/sqrt(N-3)"""
    return (1 - r**2) / np.sqrt(n - 3)

##################################################################################################################
# STEP 1: Load data
##################################################################################################################

# 1. Phenotype categories
category_df = pd.read_csv(project_DIR + '/data/unique_phenotypes_annotated.43.txt', sep='\t')
category_order = ['Metabolic', 'Hematologic', 'Hepatic', 'Musculoskeletal', 'Cardiopulmonary',  'Renal', 'Neuropsychiatric', 'Endocrine']
category_df['CATEGORY'] = pd.Categorical(category_df['CATEGORY'], categories=category_order, ordered=True)
category_df = category_df.sort_values(by=["CATEGORY", "LABEL"], ascending=[True, True])
ordered_phenotypes_abbrev = list(category_df['PHENO'])

# 2. Read couples data
# Jenny version
couples = pd.read_csv(project_DIR + '/data/hh_pairs_filter.jenny_sjaarda.txt', sep='\t')
couples = couples.astype('int32')
condition = couples['HOUSEHOLD_MEMBER1_sex'] == 1 # mask: is member1 male?
couples['HOUSEHOLD_MEMBER1'], couples['HOUSEHOLD_MEMBER2'] = np.where(condition, couples['HOUSEHOLD_MEMBER2'], couples['HOUSEHOLD_MEMBER1']), np.where(condition, couples['HOUSEHOLD_MEMBER1'], couples['HOUSEHOLD_MEMBER2'])
couples['HOUSEHOLD_MEMBER1_sex'], couples['HOUSEHOLD_MEMBER2_sex'] = 0, 1

# Stricter version
couples_strict = pd.read_csv(project_DIR+"/data/mate_pairs.yengo_def.txt", sep='\t')
couples_strict = couples_strict[['ID1', 'ID2', 'sex_id1', 'sex_id2']]
couples_strict = couples_strict.astype('int32')
condition = couples_strict['sex_id1'] == 1 # mask: is member1 male?
couples_strict['ID1'], couples_strict['ID2'] = np.where(condition, couples_strict['ID2'], couples_strict['ID1']), np.where(condition, couples_strict['ID1'], couples_strict['ID2'])
couples_strict['sex_id1'], couples_strict['sex_id2'] = 0, 1

# check no duplicate IDs
assert len(couples['HOUSEHOLD_MEMBER1']) == len(set(couples['HOUSEHOLD_MEMBER1']))
assert len(couples['HOUSEHOLD_MEMBER2']) == len(set(couples['HOUSEHOLD_MEMBER2']))
assert len(couples_strict['ID1']) == len(set(couples_strict['ID1']))
assert len(couples_strict['ID2']) == len(set(couples_strict['ID2']))

# 3. Read phenotype data: WB and from test set

# pheno_1: only test-set invidivuals 
pheno_1 = pd.read_csv(project_DIR+"/data/pheno/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.tsv.tar.gz", sep='\t')
pheno_F = pheno_1.copy()
pheno_F = pheno_F[['IID']+ordered_phenotypes_abbrev]
pheno_F.columns = [col + "_F" for col in pheno_F.columns]

pheno_M = pheno_1.copy()
pheno_M = pheno_M[['IID']+ordered_phenotypes_abbrev]
pheno_M.columns = [col + "_M" for col in pheno_M.columns]

test_IDs = set(pheno_1['IID'])

# Add phenotype data to Jenny couples
couples_v1 = couples.copy()
couples_v1 = pd.merge(couples_v1, pheno_F, left_on='HOUSEHOLD_MEMBER1', right_on='IID_F')
couples_v1 = pd.merge(couples_v1, pheno_M, left_on='HOUSEHOLD_MEMBER2', right_on='IID_M')
print("Number of couples in test + Jenny set:", len(couples_v1))

# pheno_2: All white-british individuals
pheno_2 = pd.read_csv(project_DIR+"/data/pheno/pheno_continuous_WB_INT_age_age2_sex_batch_PCs_All.tsv.tar.gz", sep='\t')
pheno_F = pheno_2.copy()
pheno_F = pheno_F[['IID']+ordered_phenotypes_abbrev]
pheno_F.columns = [col + "_F" for col in pheno_F.columns]

pheno_M = pheno_2.copy()
pheno_M = pheno_M[['IID']+ordered_phenotypes_abbrev]
pheno_M.columns = [col + "_M" for col in pheno_M.columns]

# Add phenotype data to Yengo couples
couples_v2 = couples_strict.copy()
couples_v2 = pd.merge(couples_v2, pheno_F, left_on='ID1', right_on='IID_F')
couples_v2 = pd.merge(couples_v2, pheno_M, left_on='ID2', right_on='IID_M')
print("Number of couples in WB + strict set:", len(couples_v2))

couples_v2 = couples_v2.query('ID1.isin(@test_IDs) & ID2.isin(@test_IDs)')
print("Number of couples in test + strict set:", len(couples_v2))

##################################################################################################################
# STEP 2: Compute within-couple phenotypic correlation in each set
##################################################################################################################

# Initialize empty DataFrame
corr_df = pd.DataFrame(
    index=ordered_phenotypes_abbrev,
    columns=["pheno-pheno.strict.corr",  "pheno-pheno.strict.SE", "pheno-pheno.strict.p_val", "pheno-pheno.strict.n",
             "pheno-pheno.3k.corr", "pheno-pheno.3k.SE", "pheno-pheno.3k.p_val", "pheno-pheno.3k.n"],
    dtype=float
)

# 3k
for pheno in ordered_phenotypes_abbrev:
    col_F = f"{pheno}_F"
    col_M = f"{pheno}_M"

    # Select the relevant columns and drop NA
    if col_F in couples_v1.columns and col_M in couples_v1.columns:
        valid_data = couples_v1[[col_F, col_M]].dropna()

        if not valid_data.empty:
            r, p = pearsonr(valid_data[col_F], valid_data[col_M])
            n = len(valid_data)
            se = pearson_se(r, n)

            corr_df.loc[pheno, "pheno-pheno.3k.corr"] = r
            corr_df.loc[pheno, "pheno-pheno.3k.SE"] = se
            corr_df.loc[pheno, "pheno-pheno.3k.p_val"] = p
            corr_df.loc[pheno, "pheno-pheno.3k.n"] = n

# strict
for pheno in ordered_phenotypes_abbrev:
    col_F = f"{pheno}_F"
    col_M = f"{pheno}_M"

    # Select the relevant columns and drop NA
    if col_F in couples_v2.columns and col_M in couples_v2.columns:
        valid_data = couples_v2[[col_F, col_M]].dropna()

        if not valid_data.empty:
            r, p = pearsonr(valid_data[col_F], valid_data[col_M])
            n = len(valid_data)
            se = pearson_se(r, n)

            corr_df.loc[pheno, "pheno-pheno.strict.corr"] = r
            corr_df.loc[pheno, "pheno-pheno.strict.SE"] = se
            corr_df.loc[pheno, "pheno-pheno.strict.p_val"] = p
            corr_df.loc[pheno, "pheno-pheno.strict.n"] = n

corr_df["pheno-pheno.strict.n"] = corr_df["pheno-pheno.strict.n"].astype('Int64')
corr_df["pheno-pheno.3k.n"] = corr_df["pheno-pheno.3k.n"].astype('Int64')

corr_df.to_csv('pheno_pheno.correlation.3k_vs_strict.csv')
