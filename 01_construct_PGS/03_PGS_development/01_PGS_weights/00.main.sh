# Sends all jobs for one phenotype

PHENO=$1

echo "Running pheno: $PHENO"
# 1. Run LDpred2 in parallel for different value of p (the proportion of causal variants) and capture job ID
ID_LDPRED=$(sbatch "01.run_ldpred2.sh" "$PHENO" | awk '{print $NF}')
# 2. Merge the results after LDpred2 completes
ID_MERGE=$(sbatch --dependency=afterok:"$ID_LDPRED" "02.merge_ldpred2.sh" "$PHENO" | awk '{print $NF}')


