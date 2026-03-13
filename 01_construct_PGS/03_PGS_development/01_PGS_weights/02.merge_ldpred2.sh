#!/bin/bash
#SBATCH --job-name merge_ldpred2
#SBATCH --partition urblauna
#SBATCH --output logs/%x_%A-%a.out
#SBATCH --error logs/%x_%A-%a.err
#SBATCH --partition urblauna
#SBATCH --mem 2G 
#SBATCH --time 00:10:00
#SBATCH --cpus-per-task 1
#SBATCH --get-user-env=L
#SBATCH --export NONE

#------------------------------------------------------------------------------#
# Parameters 
#------------------------------------------------------------------------------#
module load r-light

PHENO=$1
SCRATCH="/scratch/LDpred2_train/"
WORK="/data/LDpred2_train/"

# input
LDPRED_RES=$SCRATCH/LDpred2_parallel/
SUM_S_OUT=$SCRATCH/LDpred2_parallel/${PHENO}.1.pos.tsv 

# output
mkdir -p $WORK/LDpred_betas/ $WORK/LDpred_alphas/
BETAS_MERGE=$WORK/LDpred_betas/$PHENO.betas.tsv.gz
ALPHAS_MERGE=$WORK/LDpred_alphas/$PHENO.alphas.tsv

#------------------------------------------------------------------------------#
# Run Merging
#------------------------------------------------------------------------------#
Rscript Scripts/ldpred2_merge.r $PHENO $LDPRED_RES $SUM_S_OUT $ALPHAS_MERGE $BETAS_MERGE

#------------------------------------------------------------------------------#
# Add header for PGS_calc
#------------------------------------------------------------------------------#
zcat $BETAS_MERGE |\
	sed "s/rsid/rsID/g" |\
	awk 'BEGIN{print "##POLYGENIC SCORE (PGS) INFORMATION\n#genome_build=GRCh37"}{print}' |\
	gzip -c > $BETAS_MERGE.tmp
mv $BETAS_MERGE.tmp $BETAS_MERGE


