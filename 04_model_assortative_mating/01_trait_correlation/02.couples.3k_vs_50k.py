# Compare within-couple phenotypic correlation (test set vs white-british set)

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
couples = pd.read_csv(project_DIR + '/data/hh_pairs_filter.jenny_sjaarda.txt', sep='\t')
couples = couples.astype('int32')

# Condition to swap: if member1 is male, swap members so the member1 is female
condition = couples['HOUSEHOLD_MEMBER1_sex'] == 1 # mask: is member1 male?
# If condition holds, then swap members. Both operations have to be done at same time (swapping M to member2 and F to member1)
# np.where chooses first option (swap), unless it's False for condition, in which case it chooses second option (don't swap)
couples['HOUSEHOLD_MEMBER1'], couples['HOUSEHOLD_MEMBER2'] = np.where(condition, couples['HOUSEHOLD_MEMBER2'], couples['HOUSEHOLD_MEMBER1']), np.where(condition, couples['HOUSEHOLD_MEMBER1'], couples['HOUSEHOLD_MEMBER2'])
# now all first members are F, so correct sex column as well
couples['HOUSEHOLD_MEMBER1_sex'], couples['HOUSEHOLD_MEMBER2_sex'] = 0, 1

# 3. Read phenotype data: WB and test set

# pheno_1: All white-british invidivuals
pheno_1 = pd.read_csv(project_DIR+"/data/pheno/pheno_continuous_WB_INT_age_age2_sex_batch_PCs_All.tsv.tar.gz", sep='\t')
pheno_F = pheno_1.copy()
pheno_F = pheno_F[['IID']+ordered_phenotypes_abbrev]
pheno_F.columns = [col + "_F" for col in pheno_F.columns]

pheno_M = pheno_1.copy()
pheno_M = pheno_M[['IID']+ordered_phenotypes_abbrev]
pheno_M.columns = [col + "_M" for col in pheno_M.columns]

# Add phenotype data to couples df
couples_v1 = couples.copy()
couples_v1 = pd.merge(couples_v1, pheno_F, left_on='HOUSEHOLD_MEMBER1', right_on='IID_F')
couples_v1 = pd.merge(couples_v1, pheno_M, left_on='HOUSEHOLD_MEMBER2', right_on='IID_M')
print("Number of couples in WB set:", len(couples_v1))

# pheno_2: only test-set individuals
pheno_2 = pd.read_csv(project_DIR+"/data/pheno/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.tsv.tar.gz", sep='\t')
pheno_F = pheno_2.copy()
pheno_F = pheno_F[['IID']+ordered_phenotypes_abbrev]
pheno_F.columns = [col + "_F" for col in pheno_F.columns]

pheno_M = pheno_2.copy()
pheno_M = pheno_M[['IID']+ordered_phenotypes_abbrev]
pheno_M.columns = [col + "_M" for col in pheno_M.columns]

# Add phenotype data to couples df
couples_v2 = couples.copy()
couples_v2 = pd.merge(couples_v2, pheno_F, left_on='HOUSEHOLD_MEMBER1', right_on='IID_F')
couples_v2 = pd.merge(couples_v2, pheno_M, left_on='HOUSEHOLD_MEMBER2', right_on='IID_M')
print("Number of couples in test set:", len(couples_v2))

##################################################################################################################
# STEP 2: Compute within-couple phenotypic correlation in each set
##################################################################################################################

# Initialize empty DataFrame
corr_df = pd.DataFrame(
    index=ordered_phenotypes_abbrev,
    columns=["pheno-pheno.50k.corr",  "pheno-pheno.50k.SE", "pheno-pheno.50k.p_val", "pheno-pheno.50k.n",
             "pheno-pheno.3k.corr", "pheno-pheno.3k.SE", "pheno-pheno.3k.p_val", "pheno-pheno.3k.n"],
    dtype=float
)

# 50k
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

            corr_df.loc[pheno, "pheno-pheno.50k.corr"] = r
            corr_df.loc[pheno, "pheno-pheno.50k.SE"] = se
            corr_df.loc[pheno, "pheno-pheno.50k.p_val"] = p
            corr_df.loc[pheno, "pheno-pheno.50k.n"] = n

# 3k
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

            corr_df.loc[pheno, "pheno-pheno.3k.corr"] = r
            corr_df.loc[pheno, "pheno-pheno.3k.SE"] = se
            corr_df.loc[pheno, "pheno-pheno.3k.p_val"] = p
            corr_df.loc[pheno, "pheno-pheno.3k.n"] = n

corr_df["pheno-pheno.50k.n"] = corr_df["pheno-pheno.50k.n"].astype('Int64')
corr_df["pheno-pheno.3k.n"] = corr_df["pheno-pheno.3k.n"].astype('Int64')

corr_df.to_csv('pheno_pheno.correlation.50k_vs_3k.csv')

