import pandas as pd
import gzip

DIR = '/project/'

pheno_file = DIR + '/data/phenotypes/pheno_names.txt'
with open(pheno_file) as f:
    pheno_list = f.read().splitlines()

# Initialize df with IIDs only
full_df = pd.read_csv(DIR + '/data/PGS/height/score/aggregated_scores.txt.gz', sep='\t', usecols=['IID'], compression='gzip')

for pheno in pheno_list:
    print(pheno)
    PGS_df = pd.read_csv(DIR + '/data/PGS/'+pheno+'/score/aggregated_scores.txt.gz', sep='\t', compression='gzip')
    PGS_df = PGS_df[['IID', pheno+'.betas.tsv_SUM']] 
    PGS_df.rename(columns={pheno+'.betas.tsv_SUM': pheno}, inplace=True)

    # Append to full_df
    full_df = pd.merge(full_df, PGS_df, on='IID', how='left')

with gzip.open(DIR+'/data/PGS/PGS_continuous.tsv.gz', 'wt', compresslevel=9) as f:
    full_df.to_csv(f, sep='\t', index=False)

