#!/bin/bash
#SBATCH --job-name regenie_step2
#SBATCH --array=1-22
#SBATCH --output=logs/%x_%A-%a.out
#SBATCH --error=logs/%x_%A-%a.err
#SBATCH --time=3-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


CHR=${SLURM_ARRAY_TASK_ID}
threads=32

TYPE=train

mkdir -p Regenie/step2

ODIR=Regenie/step2
PFX=ukb_chr${CHR}
OUT=${ODIR}/${PFX}.${TYPE}

mkdir -p ${ODIR}

PGEN=./uk_biobank/genotypes/imp/v3/pgen/ukb22828_c${CHR}_b0_v3
COF=./UKBB/Phenotypes/Covariates.txt.gz
COL="PC{1:20},age,sex"
PRD=Regenie/step1/UKB.regenie_step1.${TYPE}_pred.list
PHF=./UKBB/Phenotypes/data/pheno_continuous_train_raw.v2.irnt.tsv.gz

# Run Regenie step 2
regenie_v3.2.9.gz_x86_64_Linux --step 2 --bsize 400 --pgen ${PGEN} --phenoFile ${PHF} --covarFile ${COF} --covarColList ${COL} --out ${OUT} --qt --ref-first --pred ${PRD} --threads ${threads} --minMAC 1

