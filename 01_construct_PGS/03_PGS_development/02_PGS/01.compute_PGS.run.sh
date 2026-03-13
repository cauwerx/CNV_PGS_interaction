#!/bin/bash

log_dir="./project/data/PGS/logs/"

# Loop through each line in the file 'pheno_list'
for pheno in $(cat  ../data/phenotypes/pheno_names.txt); do
    sbatch --output="$log_dir"/PGS_"$pheno".out \
	   --error="$log_dir"/PGS_"$pheno".err \
       00.compute_PGS.sh "$pheno"
done
