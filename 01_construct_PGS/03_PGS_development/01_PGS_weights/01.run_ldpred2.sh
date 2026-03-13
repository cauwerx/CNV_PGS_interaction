#!/bin/bash
#SBATCH --job-name LDpred2
#SBATCH --partition urblauna
#SBATCH --output logs/%x_%A-%a.out
#SBATCH --error logs/%x_%A-%a.err
#SBATCH --partition urblauna
#SBATCH --mem 60G
#SBATCH --time 2:00:00
#SBATCH --cpus-per-task 1
#SBATCH --get-user-env=L
#SBATCH --export NONE
#SBATCH --array 1-10

module load r-light

#------------------------------------------------------------------------------#
# Parameters
#------------------------------------------------------------------------------#
PHENO=$1
P_ARG=$SLURM_ARRAY_TASK_ID # proportion of causal variants. we need to test 10 different of this
SCRATCH="/scratch/LDpred2_train/"

# output
mkdir -p $SCRATCH/LDpred2_parallel/
LDPRED_RES=$SCRATCH/LDpred2_parallel/${PHENO}.${P_ARG}.multi_auto.rds 
SUM_S_OUT=$SCRATCH/LDpred2_parallel/${PHENO}.${P_ARG}.pos.tsv # save the filtered positions

#------------------------------------------------------------------------------#
# Run
#------------------------------------------------------------------------------#
Rscript Scripts/ldpred2.r $PHENO $P_ARG $LDPRED_RES $SUM_S_OUT

