#!/bin/bash
#SBATCH --job-name regenie
#SBATCH --array=0-0
#SBATCH --output=logs/%x_%A-%a.out
#SBATCH --error=logs/%x_%A-%a.err
#SBATCH --time=3-0
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

TYPE=train

mkdir -p Regenie/step1

DDIR=./UKBB/GWAS/data
PDIR=./UKBB/Phenotypes
BED=${DDIR}/_001_ukb_cal
SNP=${DDIR}/qc_pass.snplist
SMP=${DDIR}/qc_pass.id
PHF=${PDIR}/data/pheno_continuous_train_raw.v2.irnt.tsv.gz
COF=${PDIR}/Covariates.txt.gz
COL="PC{1:20},age,sex"
OUT=Regenie/step1/UKB.regenie_step1.${TYPE}

threads=32

# Run Regenie step 1
regenie_v3.2.9.gz_x86_64_Linux --step 1 --bsize 1000 --bed ${BED} --extract ${SNP} --keep ${SMP} --phenoFile ${PHF} --covarFile ${COF} --covarColList ${COL} --lowmem --lowmem-prefix Regenie/step1/regenie_tmp_preds.${TYPE} --qt --out ${OUT} --threads ${threads}

