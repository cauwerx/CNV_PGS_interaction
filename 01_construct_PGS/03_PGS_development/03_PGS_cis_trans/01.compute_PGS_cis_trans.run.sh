#!/bin/bash

log_dir="./project/data/PGS.in_out/logs/"

# Loop through each line in the file 'pheno_list'
for pheno in $(cat  ../data/phenotypes/pheno_names.with_CNV_signal.txt); do
    sbatch --output="$log_dir"/PGS_"$pheno".out \
	   --error="$log_dir"/PGS_"$pheno".err \
       01.compute_PGS_cis_trans.sh "$pheno"
done


