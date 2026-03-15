# Save info on phenotype (with or without location/sex correction), PGS_GW or PGS_trans and CNV status

################################################################################
# Libraries
################################################################################
import pandas as pd
import numpy as np

project_DIR = "/mnt/project/"

################################################################################
# STEP 1: Load phenotype, PGS and CNV info
################################################################################

# 1.1) INT+cov-corrected phenotypes and PGSs
pheno = pd.read_csv(project_DIR + "/data/pheno/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.tsv.tar.gz", sep='\t')
#pheno = pd.read_csv(project_DIR + "/data/pheno_with_location_corr/pheno_continuous_test_INT_age_age2_sex_batch_location_PCs_All.tsv.tar.gz", sep='\t') # with coordinates correction
#pheno = pd.read_csv(project_DIR + "/data/pheno_with_sex_interaction_correction/pheno_continuous_test_INT_age_age2_sex_batch_PCs_All.tsv", sep='\t') # with sex interaction correction
pgs = pd.read_csv(project_DIR + "/data/PGS.out_of_sample/PGS_continuous.tsv.gz", sep='\t') # PGS_GW
#pgs = pd.read_csv(project_DIR + "/data/PGS.out_of_sample/PGS_continuous.out.250kb.tsv.gz", sep='\t') # PGS_trans

print('Number of individuals with phenotype:', len(pheno))
print('Number of individuals with PGS:', len(pgs))

# Change PGS column names to distinguish them from phenotype columns
pgs.columns = ['IID'] + [col + '_PGS' for col in pgs.columns.drop('IID')]

# Merge pheno with PGS 
full_info = pheno.copy()
full_info = full_info.merge(pgs, on="IID", how="left")

print('Number of individuals with PGS & pheno:', len(full_info))

# 1.2) CNV carriers
cnv = pd.read_table(project_DIR + "/data/ukb_cnv_global.gz")
# Obtain IID and convert to int
cnv['IID'] = cnv['Sample_Name'].apply(lambda x: x.split('_')[0])
cnv['IID'] = cnv['IID'].astype(int)
cnv['Chromosome'] = cnv['Chromosome'].astype(int)
cnv['Start_Position_bp'] = cnv['Start_Position_bp'].astype(int)
cnv['End_Position_bp'] = cnv['End_Position_bp'].astype(int)

print('Number of individuals with CNV calls:', len(set(cnv['IID'])))

# CNV signals (CNVR associated to a complex trait) 
CNVR = pd.read_table(project_DIR + "/data/CNV_GWAS_signals_ranking_v1.txt", sep="\t")
# Select only continuous phenotypes, not sex-specific
CNVR = CNVR.query('TYPE=="continuous" & SEX=="All"')
# Remove chromosome X
CNVR = CNVR[CNVR['CHR']!='X']
CNVR['REGION'] = CNVR['CHR']+'_'+CNVR['CNVR_START'].astype(str)+'_'+CNVR['CNVR_STOP'].astype(str)
CNVR['CHR'] = CNVR['CHR'].astype(int)
CNVR['TOP_POS'] = CNVR['TOP_POS'].astype(int)
CNVR = CNVR.sort_values(by=['PHENO', 'CHR', 'CNVR_START'])
# There is one top model which is both mirror (M) and duplication (D). Use M since it's more general.
both_model_i = CNVR.query('TOP_MODEL=="M-DUP"').index[0]
CNVR.loc[both_model_i, 'TOP_MODEL']="M"

print('Number of CNV signals:', CNVR.shape[0])

################################################################################
# STEP 2: Identify carriers of trait-associated CNVs and encode them
################################################################################

# CNV carriers defined as containing to lead probe (with location=chrom:position)
def get_CNV_carriers(chrom, position, cnv): # cnv is ukb_cnv_global.gz file
    # High-quality (abs(QS)>0.5) 
    # DEL carriers
    del_carriers=set(cnv.query('Quality_Score<-0.5 & Chromosome==@chrom & Start_Position_bp<=@position & End_Position_bp>=@position')['IID'])
    # DUP carriers
    dup_carriers=set(cnv.query('Quality_Score>0.5 & Chromosome==@chrom & Start_Position_bp<=@position & End_Position_bp>=@position')['IID'])
    # Low quality carriers
    lq_carriers=set(cnv.query('abs(Quality_Score)<0.5 & Chromosome==@chrom & Start_Position_bp<=@position & End_Position_bp>=@position')['IID'])
    
    return(del_carriers, dup_carriers, lq_carriers)

# Encode according to association model and CNV status:
# Mirror model: (DEL, CN, DUP) = (-1, 0, 1)
# DEL model: (DEL, CN) = (1, 0)
# DUP model: (CN, DUP) = (0, 1)

# Initialize storing df
CNV_status = pheno[['IID']].copy()

# Iterate over each pheno-CNVR pair
for index, row in CNVR.iterrows():    
    phenotype = row['PHENO']
    chrom = row['CHR']
    pos = row['TOP_POS']
    region = row['REGION']
    model = row['TOP_MODEL']
    CNV_name = phenotype + '_CNV_' + region # these values are unique
    print("Processing: ", CNV_name)
    
    # Get carriers in that region (containing the lead probe)
    # Notice this set contains all UKBB participants. In the pheno df, we only have a subset of them (test set)
    del_carriers, dup_carriers, lq_carriers = get_CNV_carriers(chrom, pos, cnv)
    
    # By default set carrier status to zero in that CNVR for that phenotype. Each CNV_name is a different column.
    CNV_status[CNV_name]=0
    
    if model=="DEL":
        CNV_status.loc[CNV_status['IID'].isin(del_carriers), CNV_name]=1
        CNV_status.loc[CNV_status['IID'].isin(dup_carriers), CNV_name]=np.nan
        CNV_status.loc[CNV_status['IID'].isin(lq_carriers),  CNV_name]=np.nan
        # Check DEL carriers in full set are more than in test subset
        assert(len(CNV_status.loc[CNV_status[CNV_name]==1]) <= len(del_carriers))
    
    elif model=="DUP":
        CNV_status.loc[CNV_status['IID'].isin(del_carriers), CNV_name]=np.nan
        CNV_status.loc[CNV_status['IID'].isin(dup_carriers), CNV_name]=1
        CNV_status.loc[CNV_status['IID'].isin(lq_carriers),  CNV_name]=np.nan 
        # Check DUP carriers in full set are more than in test subset
        assert(len(CNV_status.loc[CNV_status[CNV_name]==1]) <= len(dup_carriers))
               
    elif model=="M":
        CNV_status.loc[CNV_status['IID'].isin(del_carriers), CNV_name]=-1
        CNV_status.loc[CNV_status['IID'].isin(dup_carriers), CNV_name]=1
        CNV_status.loc[CNV_status['IID'].isin(lq_carriers),  CNV_name]=np.nan
        # Check DEL carriers in full set are more than in test subset
        assert(len(CNV_status.loc[CNV_status[CNV_name]==-1]) <= len(del_carriers))
        # Check DUP carriers in full set are more than in test subset
        assert(len(CNV_status.loc[CNV_status[CNV_name]==1]) <= len(dup_carriers))
               
    else:
        print('Model not available.')

        
# Some checks
cnv_cols = CNV_status.loc[:, CNV_status.columns.str.contains("CNV_")].columns
# Check columns are unique
assert len(cnv_cols)==len(set(cnv_cols))
# Check that number of CNV columns = number of CNV-trait pairs
assert len(cnv_cols)==len(CNVR)

################################################################################
# STEP 3: Merge CNV_status with pheno and PGS
################################################################################

# Check the full_info df and CNV_status df have same IID number
assert len(full_info)==len(CNV_status)
full_info = full_info.merge(CNV_status, on="IID", how="left")

################################################################################
# STEP 4: Save pheno + PGS + CNV_status info
################################################################################

full_info.to_csv('pheno_pgs_cnv.autosomes.non_sex_specific.csv', index=False) # PGS_GW
#full_info.to_csv('pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.csv', index=False) # PGS_trans
#full_info.to_csv('pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.pheno_location_correction.csv', index=False) # pheno with coordinates correction
#full_info.to_csv('pheno_pgs_cnv.autosomes.non_sex_specific.pgs_trans.pheno_sex_correction.csv', index=False) # pheno with sex correction
