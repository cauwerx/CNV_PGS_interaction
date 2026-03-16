#!/bin/bash
#SBATCH --job-name bgen2pgen
#SBATCH --output /data/logs/pgen_chr_%a.out
#SBATCH --error  /data/logs/pgen_chr_%a.err
#SBATCH --partition urblauna
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 500G
#SBATCH --time 03:30:00
#SBATCH --get-user-env=L
#SBATCH --export NONE
#SBATCH --array=1-22

module load plink2/2.00a4.3

CHR=$SLURM_ARRAY_TASK_ID
BGEN=/data/uk_biobank/genotypes/imp/v3/ukb22828_c${CHR}_b0_v3.bgen
SAMPLES=/data/uk_biobank/genotypes/imp/ukb1638_imp_chr1_v2_s487398.sample
PGEN=/scratch/pgen_v3/ukb22828_c${CHR}_b0_v3

plink2 --bgen $BGEN ref-first --sample $SAMPLES --mach-r2-filter 0.3 2 --maf 0.005 --make-pgen --out $PGEN
