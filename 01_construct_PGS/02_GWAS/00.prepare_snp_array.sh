#!/bin/bash
#SBATCH --job-name merge_array
#SBATCH --array=0-0
#SBATCH --output=logs/%x_%A-%a.out
#SBATCH --error=logs/%x_%A-%a.err
#SBATCH --time=1-0
#SBATCH --cpus-per-task=48
#SBATCH --mem=500gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

rm -f list_beds.txt
touch list_beds.txt

# 1. Location of PLINK files split by CHR number
BDIR=./plinkfiles_DIR
B_PFX=_001_ukb_cal

# 2. Create list of file paths (bed, bim, fam) from CHR 2 to 22
for CHR in {2..22}; do
	PFX=${B_PFX}_chr${CHR}_v2
	echo "${BDIR}/${PFX}.bed ${BDIR}/${PFX}.bim ${BDIR}/${PFX}.fam" >> list_beds.txt
done

# 3. Set CHR 1 as basefile for merging
CHR=1
PFX=${B_PFX}_chr${CHR}_v2
BED=${BDIR}/${PFX}.bed
BIM=${BDIR}/${PFX}.bim
FAM=${BDIR}/${PFX}.fam
OUT=Out/SNParray/${B_PFX}
THR=48

mkdir -p data

# 4. Merge all autosomes with PLINK1
plink --bed ${BED} --bim ${BIM} --fam ${FAM} --merge-list list_beds.txt --make-bed --out ${OUT} --threads ${THR}

# 5. Perform variant and sample QC with PLINK2
plink2 \
  --bfile ${OUT} \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --write-snplist --write-samples --no-id-header \
  --out data/qc_pass


