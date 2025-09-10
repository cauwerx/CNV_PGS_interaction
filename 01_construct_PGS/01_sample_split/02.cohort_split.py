# Split white-british set into non-carrier training and carrier-enriched test splits

################################################################################
# Libraries
################################################################################

# Install libraries
#!pip install -U scikit-learn
#!pip install pyarrow

# Load libraries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import pyarrow

HOME_DIR = "/mnt/project/penetrance_at_scale/"

################################################################################
# STEP 1: Load full cohort and select white-british participants
################################################################################

# sample_eid: this corresponds to a list of all sample eids (with sex information), retrieved from a ".fam" file
# 0=eid, 4=sex
sample_eid = pd.read_csv("/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c1_b0_v2.fam", sep=' ', usecols=[0,4], header = None, names=["eid", "sex"])

# sample_QC: this corresponds to the "ukb_sqc_v2.txt" file 
# Get column names obtained in 0_sampleQC_colnames.R script
f = open(HOME_DIR + "/data/helper_files/ukb_sqc_v2.HEADER.txt")
col_names = f.read().splitlines() 
# Load the QC file
sample_qc = pd.read_csv("/mnt/project/Bulk/Genotype Results/Genotype calls/ukb_sqc_v2.txt", sep=' ', header=None, names=col_names)

# Merge both column-wise
full_cohort = pd.concat([sample_eid, sample_qc], axis=1)

# Select white-british individuals
WB_cohort = full_cohort.query('in_white_british_ancestry_subset==1').copy()
WB_IDs = set(WB_cohort['eid'])

# Check number of WB participants
print("Numer of white-british individuals:", len(WB_IDs))

################################################################################
# STEP 2: Load CNV calls, CNV-burden calls and CNV-trait associations
################################################################################

# Import CNV carriers
cnv = pd.read_table(HOME_DIR + "/data/ukb_cnv_global.gz")
# Obtain IID and convert to int
cnv['IID'] = cnv['Sample_Name'].apply(lambda x: x.split('_')[0])
cnv['IID'] = cnv['IID'].astype(int)
cnv['Chromosome'] = cnv['Chromosome'].astype(int)
cnv['Start_Position_bp'] = cnv['Start_Position_bp'].astype(int)
cnv['End_Position_bp'] = cnv['End_Position_bp'].astype(int)

# Import CNV-burden
dtype_dict = {
    "IID": int,
    "SEX": "category",
    "BURDEN_QS": float,
    "BURDEN_Mb": float,
    "BURDEN_Num100": float,
    "BURDEN_Num250": float,
    "BURDEN_Num500": float,
    "BURDEN_Num1000": float,
    "BURDEN_Num2000": float,
    "BURDEN_GENES": float
}
cnv_burden=pd.read_csv(HOME_DIR +'/data/CNV_burden/CNV_burden.txt.gz', compression='gzip', sep='\t', dtype=dtype_dict)

# Import CNV signals (CNVR associated to a phenotype)
CNVR = pd.read_table(HOME_DIR + "data/CNV_GWAS_signals_ranking_v1.txt", sep="\t")
CNVR = CNVR.sort_values(by='PHENO')
CNVR['CHR'] = CNVR['CHR'].replace({'X': '23'})
CNVR['CHR'] = CNVR['CHR'].astype(int)
CNVR['CNVR_START'] = CNVR['CNVR_START'].astype(int)
CNVR['CNVR_STOP'] = CNVR['CNVR_STOP'].astype(int)
# There is one top model which is both DUP and M. Select M model (more general).
twoM_index = CNVR.query('TOP_MODEL=="M-DUP"').index[0]
CNVR.loc[twoM_index,'TOP_MODEL'] = "M"

################################################################################
# STEP 3: Select high-quality CNV carriers & high-quality CNV-burden carriers
################################################################################

# 3.1) High-quality CNV carriers
# CNV carriers defined as containing to lead probe
def get_CNV_carriers(chrom, position, cnv): #cnv is ukb_cnv_global.gz
    # High-quality (abs(QS)>0.5) 
    # DEL carriers
    hq_del_carriers=list(set(cnv.query('Quality_Score<-0.5 & Chromosome==@chrom & Start_Position_bp<=@position & End_Position_bp>=@position')['IID']))
    # DUP carriers
    hq_dup_carriers=list(set(cnv.query('Quality_Score> 0.5 & Chromosome==@chrom & Start_Position_bp<=@position & End_Position_bp>=@position')['IID']))
    
    # Low-quality (0.2<abs(QS)<=0.5) 
    # low quality DEL carriers
    lq_del_carriers=list(set(cnv.query('Quality_Score<-0.2 & Quality_Score>=-0.5 & Chromosome==@chrom & Start_Position_bp<=@position & End_Position_bp>=@position')['IID']))
    # low quality DUP carriers
    lq_dup_carriers=list(set(cnv.query('Quality_Score> 0.2 & Quality_Score<= 0.5 & Chromosome==@chrom & Start_Position_bp<=@position & End_Position_bp>=@position')['IID']))
    
    return(hq_del_carriers, hq_dup_carriers, lq_del_carriers, lq_dup_carriers)

# Save high-quality phenotype-assocaited CNV carriers
CNV_carriers = []
LQ_CNV_carriers = []
for i, row in CNVR.iterrows():
    chrom=row['CHR']
    position=row['TOP_POS']
    model=row['TOP_MODEL']
    hq_del_carriers, hq_dup_carriers, lq_del_carriers, lq_dup_carriers = get_CNV_carriers(chrom, position, cnv)
    if model=='DEL':
        CNV_carriers.extend(hq_del_carriers)
        LQ_CNV_carriers.extend(lq_del_carriers)
    elif model=='DUP':
        CNV_carriers.extend(hq_dup_carriers)
        LQ_CNV_carriers.extend(lq_dup_carriers)
    elif model=='M':
        CNV_carriers.extend(hq_del_carriers)
        CNV_carriers.extend(hq_dup_carriers)
        LQ_CNV_carriers.extend(lq_del_carriers)
        LQ_CNV_carriers.extend(lq_dup_carriers)

# Remove repeated IDs
CNV_carriers = set(CNV_carriers)
LQ_CNV_carriers = set(LQ_CNV_carriers)

# 3.2) High-quality CNV-burden carriers (genes_affected >=10)
high_burden_carriers = list(cnv_burden.query('BURDEN_GENES>=10')['IID'].values)

print("Number of high quality carriers:", len(CNV_carriers))
print("Number of low quality carriers:", len(LQ_CNV_carriers))
print("Number of high burden carriers:", len(high_burden_carriers))

# 3.3) Join the CNV_carriers set with high CNV burden carriers
all_carriers = CNV_carriers.union(high_burden_carriers)

print("Number of HQ carriers + high-burden carriers:", len(CNV_carriers)+len(high_burden_carriers))
# Some individuals are shared between the sets
print("Number of all carriers:", len(all_carriers))

# 3.4) Keep the carriers that are white-british
all_carriers = all_carriers.intersection(WB_IDs)
print("Number of WB all carriers:", len(all_carriers))

################################################################################
# STEP 4: Remove relatives of carriers
################################################################################

# Import relatedness info
dtype_dict = {
    "ID1": int,
    "ID2": int,
    "HetHet": float,
    "IBS0": float,
    "Kinship": float
}
relatedness = pd.read_csv("/mnt/project/Bulk/Genotype Results/Genotype calls/ukb_rel.dat", sep=' ', dtype=dtype_dict)

# Get at least 3rd degree relatives (or more related) that are WB
related = relatedness.query('Kinship>0.0884 & ID1.isin(@WB_IDs) & ID2.isin(@WB_IDs)').copy()


# Remove carriers with the highest numbers of relatives in the test set
test_exclusion_lst = []

# Get the relative pairs where both are carriers
tmp_df = related.query('ID1.isin(@all_carriers) & ID2.isin(@all_carriers)').copy()

# While I haven't finished reading the tmp_df (which shrinks progressively)...
while len(tmp_df) != 0:
    # Find person with the most relatives:
    # Combine both ID columns into a single Series
    all_ids = pd.concat([tmp_df['ID1'], tmp_df['ID2']])
    # Count how many times each ID appears, the index is now the ID
    id_counts = all_ids.value_counts()
    # Get the person (index in Series) who appears the most (highest value)
    most_common_person = id_counts.idxmax()

    # If person isn't a carrier raise an error (as both should be carriers)
    if most_common_person not in all_carriers:
        print("Non carrier found")
        break
    # Save person and remove person from tmp_df
    test_exclusion_lst.append(most_common_person)
    tmp_df = tmp_df.query('ID1!=@most_common_person & ID2!=@most_common_person').copy()

test_exclusion_lst = set(test_exclusion_lst)

# Save all relatives in all_carriers
related_all_carriers = all_carriers

# Remove relatives among all_carriers
all_carriers = all_carriers.difference(test_exclusion_lst)

print("Number of non-related WB all carriers:", len(all_carriers))

################################################################################
# STEP 5: Removing carriers and their relatives from training set
################################################################################

# Select white-british individuals
train_IDs = WB_IDs

# Remove all HQ CNV carriers and high burden individuals, as well as LQ CNV carriers
train_IDs = train_IDs.difference(related_all_carriers) # we want to rm all carriers, not only non-related ones
train_IDs = train_IDs.difference(LQ_CNV_carriers)
print("Number of non-carrier individuals:", len(train_IDs))

# Find non-carriers related to "all_carriers" set and remove them from training set
# Only one of the pair is in "all_carriers"
tmp_df = related.query('(ID1 in @all_carriers & ID2 in @train_IDs) | (ID2 in @all_carriers & ID1 in @train_IDs)').copy()
# Who is in "all_carriers" set?
mask = tmp_df['ID1'].isin(all_carriers)
# Swap them such that ID1 is always the carrier (this operation needs to be done simultaneously!)
tmp_df['ID1'], tmp_df['ID2'] = np.where(mask, tmp_df['ID1'], tmp_df['ID2']), np.where(mask, tmp_df['ID2'], tmp_df['ID1'])
# Save exclusion set (non-carriers)
train_exclusion_set = set(tmp_df['ID2'])

# Remove non-carrier relatives of "all_carriers"
train_IDs = train_IDs.difference(train_exclusion_set)
print("Number of non-carrier individuals, not related to carriers:", len(train_IDs))

################################################################################
# STEP 6: Dividing non-carrier set into training set and control set
################################################################################

# Load birth_year info
age_df = pd.read_parquet(HOME_DIR + '/data/pheno/age.all.parquet', engine='pyarrow')
# Get birth_year and sex for WB individuals
train_merged = pd.merge(age_df[['IID','birth_year']], WB_cohort[['eid','sex']], how='inner', left_on='IID', right_on='eid')
train_merged.drop(columns='eid', inplace=True)
# Filter for individuals in training set
train_merged = train_merged.query('IID.isin(@train_IDs)').copy()
# Convert birth_year to age (assuming current year = 2025)
train_merged['age'] = 2025 - train_merged['birth_year']

# Bin by age
train_merged['age_quintile'] = pd.qcut(train_merged['age'], q=5, labels=False) # discretize based on age quintiles, labels=False to return integer bins
# Add sex info to age bins
train_merged['strata'] = train_merged['sex'].astype(str) + '_' + train_merged['age_quintile'].astype(str)
# Stratified split based on age_sex bin (divides into training and test with same proportions for "strata")
train_df, test_df = train_test_split(train_merged, test_size=0.2, stratify=train_merged['strata'], random_state=42)

# Actual training set IDs
train_IDs = set(train_df['IID'])
print("Number of individuals in training set:", len(train_IDs))

# Test set IDs: carriers + some "control" non_carriers
test_IDs = all_carriers.union(set(test_df['IID']))
print("Number of individuals in test set:", len(test_IDs))

# Perform some checks
# Sex percentage:
summary = train_df['sex'].value_counts(normalize=True).mul(100).round(1).astype(str) + '%'
print("Train")
print(summary)
summary = test_df['sex'].value_counts(normalize=True).mul(100).round(1).astype(str) + '%'
print("Control Test")
print(summary)
print('\n')

# Age quintile percentage:
summary = train_df['age_quintile'].value_counts(normalize=True).mul(100).round(1).astype(str) + '%'
print("Train")
print(summary)
summary = test_df['age_quintile'].value_counts(normalize=True).mul(100).round(1).astype(str) + '%'
print("Control Test")
print(summary)

################################################################################
# STEP 7: Remove relatives among non-carriers in test set
################################################################################

# Individuals in all_carriers are not related
assert related.query('ID1.isin(@all_carriers) & ID2.isin(@all_carriers)').shape[0] == 0, 'Relatives found in all_carriers'

# Individuals in all_carriers are not related to control non-carriers
assert related.query('ID1.isin(@all_carriers) & ID2.isin(@test_IDs)').shape[0] == 0, 'Relatives between all_carriers and control non-carriers'
assert related.query('ID1.isin(@test_IDs) & ID2.isin(@all_carriers)').shape[0] == 0, 'Relatives between all_carriers and control non-carriers'

# There are a few non-carriers that still have other non-carriers relatives in the test set
# Remove non-carriers with the highest numbers of relatives in the test set
test_exclusion_lst=[]

# Get the relatives pairs where both are non-carriers and both belong to the test set
tmp_df = related.query('ID1 in @test_IDs & ID2 in @test_IDs & ID1 not in @all_carriers & ID2 not in @all_carriers').copy()

# While I haven't finished reading the tmp_df (which shrinks progressively)...
while len(tmp_df) != 0:
    # Find person with the most relatives:
    # Combine both ID columns into a single Series
    all_ids = pd.concat([tmp_df['ID1'], tmp_df['ID2']])
    # Count how many times each ID appears, the index is now the ID
    id_counts = all_ids.value_counts()
    # Get the person (index in Series) who appears the most (highest value)
    most_common_person = id_counts.idxmax()

    # If person is a carrier raise an error (as both shouldn't be carriers)
    if most_common_person in all_carriers:
        print("Carrier found")
        break
    # Save person and remove person from tmp_df
    test_exclusion_lst.append(most_common_person)
    tmp_df = tmp_df.query('ID1!=@most_common_person & ID2!=@most_common_person').copy()

test_exclusion_lst = set(test_exclusion_lst)

# Update test IDs
test_IDs = test_IDs.difference(test_exclusion_lst)

print("Number of non-related individuals in test set:", len(test_IDs))

# Check all relatives were removed
assert related.query('ID1.isin(@test_IDs) & ID2.isin(@test_IDs)').shape[0] == 0, 'Relatives found in test_IDs'

################################################################################
# STEP 8: Remove relatives between non-carriers in training and test sets
################################################################################

# There are still some relatives between training and test set
# Find non-carriers in test set related to those non-carriers in training set
# One of the pair is in test set and the other in training set
tmp_df = related.query('(ID1 in @test_IDs & ID2 in @train_IDs) | (ID2 in @test_IDs & ID1 in @train_IDs)').copy()
# Who is the one in test set?
mask = tmp_df['ID1'].isin(test_IDs)
# Swap them such that ID1 is always the one in test (this operation needs to be done simultaneously!)
tmp_df['ID1'], tmp_df['ID2'] = np.where(mask, tmp_df['ID1'], tmp_df['ID2']), np.where(mask, tmp_df['ID2'], tmp_df['ID1'])
# Save the IDs in test
test_exclusion_lst = set(tmp_df['ID1'])
# Remove these IDs from test
test_IDs = test_IDs.difference(test_exclusion_lst)

print("Number of non-related individuals in test set, also not related to individuals in training set:", len(test_IDs))

# Check all relatives were removed
assert related.query('ID1.isin(@test_IDs) & ID2.isin(@train_IDs)').shape[0] == 0, 'Relatives found between train_IDs & test_IDs'
assert related.query('ID1.isin(@train_IDs) & ID2.isin(@test_IDs)').shape[0] == 0, 'Relatives found between train_IDs & test_IDs'

################################################################################
# STEP 9: Save IDs
################################################################################

trainIDs_df = pd.DataFrame(train_IDs, columns=["IID"])
testIDs_df = pd.DataFrame(test_IDs, columns=["IID"])

trainIDs_df.to_csv('train_IDs.txt', index=False, header=True, sep='\t')
testIDs_df.to_csv('test_IDs.txt', index=False, header=True, sep='\t')

