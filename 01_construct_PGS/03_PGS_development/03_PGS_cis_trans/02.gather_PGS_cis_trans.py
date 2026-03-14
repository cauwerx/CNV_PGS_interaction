################################################################################
# Libraries
################################################################################
import pandas as pd
import gzip

################################################################################
# Parameters
################################################################################

DIR = './project'
pheno_file = DIR + '/data/phenotypes/pheno_names.with_CNV_signal.txt'
with open(pheno_file) as f:
    pheno_list = f.read().splitlines()

################################################################################
# STEP 1: Get IIDs and select window type
################################################################################

# Initialize df (only IID)
full_in_df = pd.read_csv(DIR + '/data/PGS.in_out/height/score/aggregated_scores.txt.gz', sep='\t', usecols=['IID'], compression='gzip')
full_out_df = pd.read_csv(DIR + '/data/PGS.in_out/height/score/aggregated_scores.txt.gz', sep='\t', usecols=['IID'], compression='gzip')

size = "" #50 kb, 100kb, 250kb, empty string for CHR

################################################################################
# STEP 2: Merge all same window-size files across phenotypes
################################################################################

# in = cis, out = trans

for pheno in pheno_list:
    print(pheno)
    PGS_df = pd.read_csv(DIR + '/data/PGS.in_out/'+pheno+'/score/aggregated_scores.txt.gz', sep='\t', compression='gzip')

    if 'kb' in size:
        # PGS_in, select and rename cols
        PGSin_df = PGS_df[['IID'] + [col for col in PGS_df.columns if col.endswith('_'+size+'_in_SUM')]].copy()
        PGSin_df.rename(columns={col: col.replace('_'+size+'_in_SUM', '') for col in PGSin_df.columns if col.endswith('_'+size+'_in_SUM')}, inplace=True)

        # PGS_out, select and rename cols
        PGSout_df = PGS_df[['IID'] + [col for col in PGS_df.columns if col.endswith('_'+size+'_out_SUM')]].copy()
        PGSout_df.rename(columns={col: col.replace('_'+size+'_out_SUM', '') for col in PGSout_df.columns if col.endswith('_'+size+'_out_SUM')}, inplace=True)

    elif size=="":
        # PGS_in, select and rename cols
        PGSin_df = PGS_df[['IID'] + [col for col in PGS_df.columns if col.endswith('_in_SUM') and 'kb' not in col]].copy()
        PGSin_df.rename(columns={col: col.replace('_in_SUM', '') for col in PGSin_df.columns if col.endswith('_in_SUM')}, inplace=True)

        # PGS_out, select and rename cols
        PGSout_df = PGS_df[['IID'] + [col for col in PGS_df.columns if col.endswith('_out_SUM') and 'kb' not in col]].copy()
        PGSout_df.rename(columns={col: col.replace('_out_SUM', '') for col in PGSout_df.columns if col.endswith('_out_SUM')}, inplace=True)

    else:
        print("Not available.")
        break

    # Append to full_df (for all phenotypes)
    full_in_df = pd.merge(full_in_df, PGSin_df, on='IID', how='left')
    full_out_df = pd.merge(full_out_df, PGSout_df, on='IID', how='left')

if size == "":
    window = "CHR"
else:
    window = size

################################################################################
# STEP 3: Save cis and trans PGS for given window-size separately
################################################################################

with gzip.open(DIR+'/data/PGS.in_out/PGS_continuous.in.'+window+'.tsv.gz', 'wt', compresslevel=9) as f_in:
    full_in_df.to_csv(f_in, sep='\t', index=False)

with gzip.open(DIR+'/data/PGS.in_out/PGS_continuous.out.'+window+'.tsv.gz', 'wt', compresslevel=9) as f_out:
    full_out_df.to_csv(f_out, sep='\t', index=False)

