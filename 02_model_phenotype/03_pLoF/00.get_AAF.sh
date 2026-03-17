path_to_500kwes_plink_files="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/"

for CHR in {1..22} X Y; do
    dx run app-swiss-army-knife \
        -iin="${path_to_500kwes_plink_files}/ukb23158_c${CHR}_b0_v1.bed" \
        -iin="${path_to_500kwes_plink_files}/ukb23158_c${CHR}_b0_v1.bim" \
        -iin="${path_to_500kwes_plink_files}/ukb23158_c${CHR}_b0_v1.fam" \
        -iin="/project/data/samples_white_british_All.txt" \
        -icmd="plink2 --bfile ukb23158_c${CHR}_b0_v1 --no-fam-pheno --keep samples_white_british_All.txt --freq --out ukb23158_c${CHR}_b0_v1.nonrelated_whitebritish" \
        --destination "/project/data/helper_files/aaf_files/" \
        --instance-type mem1_ssd1_v2_x16 --priority normal \
        --name aaf_chr${CHR}_all_variants_wb --tag aaf_all_w -y
done
