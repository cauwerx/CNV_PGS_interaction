import pandas as pd

DIR='/project/data/'

pheno_map = pd.read_csv(DIR+'/phenotypes/pheno_list.txt', sep='\t')

# Load list of final phenotypes (some pheno where computed rather than extracted)
actualpheno_file = DIR+'/phenotypes/pheno_continuous_WB_raw.tsv'
with open(actualpheno_file, 'r') as f:
    actualpheno_list = f.readline().strip().split('\t')

# Remove non-pheno items from list
not_pheno_list = ['IID', 'sex']
actualpheno_list = [item for item in actualpheno_list if item not in not_pheno_list]
# Pheno that were computed don't have a FieldID
no_ID_list = ['GS', 'WHR', 'WHRadjBMI']

for pheno_name in actualpheno_list:
    if pheno_name not in no_ID_list:
        pheno_ID = pheno_map.query('Pheno==@pheno_name')['FieldID'].iloc[0]
    else: 
        pheno_ID = '-'
    print(pheno_name, pheno_ID)

    score_df = pd.read_csv(DIR+'/PGS_scorefiles/' + pheno_name + '.betas.tsv.gz', sep='\t', skiprows=2, compression='gzip')

    # Re-arrange columns
    score_df = score_df[['chr_name', 'chr_position', 'rsID', 'effect_allele', 'other_allele', 'effect_weight']]
    
    newscore_file = DIR+'/PGS_scorefiles.parsed/' + pheno_name + '.betas.tsv'

    with open(newscore_file, 'w') as f:
        f.write(f"#pgs_name={pheno_name}\n")
        f.write(f"#pgs_id={pheno_ID}\n")
        f.write(f"#trait_reported={pheno_name}\n")
        f.write(f"#genome_build=GRCh37\n")

    score_df.to_csv(newscore_file, sep='\t', index=False, mode='a')

