#!/bin/bash
#SBATCH --job-name getPGS_in_out
#SBATCH --partition urblauna
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 100G
#SBATCH --time 40:00:00
#SBATCH --get-user-env=L
#SBATCH --export=NONE

# Get phenotype from command line argument, default to "height"
pheno=${1:-height}

# Run pgsc_calc
module load plink2/2.00a4.3
DIR=./project/

${DIR}/software/nextflow-23.04.2-all run ${DIR}/software/pgsc_calc/main.nf \
  -profile singularity \
  -work-dir ${DIR}/data/PGS.in_out/${pheno}/work/ \
  --input ${DIR}/data/pgen_path.v3.csv \
  --scorefile "${DIR}/data/PGS_scorefiles.parsed.in_out/${pheno}_chr*.txt" \
  --target_build GRCh37 \
  --ref ${DIR}/software/pgsc_calc_ref.sqlar \
  --outdir ${DIR}/data/PGS.in_out/${pheno}/

