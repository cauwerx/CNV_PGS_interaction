# Compute pLoF burden with PLINK

type="pLoF"
path_to_score_files="/penetrance_at_scale/data/scorefiles_perCHR_MAF0.01/"
path_to_500kwes_plink_files="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/"

for scorefile in $(cat ./data/scorefiles_perCHR_MAF0.01/scorefiles.${type}.list)
do 
  basename=$(basename "$scorefile" .txt)
  chr_num=$(echo "$basename" | cut -d_ -f1 | grep -oE '[0-9]+')
  if [[ $chr_num -eq 23 ]]; then
    chr_num="X"
  elif [[ $chr_num -eq 24 ]]; then
    chr_num="Y"
  fi
  dx run swiss-army-knife \
    -iin="${path_to_score_files}/${basename}.txt" \
  	-iin="${path_to_500kwes_plink_files}/ukb23158_c${chr_num}_b0_v1.bed" \
    -iin="${path_to_500kwes_plink_files}/ukb23158_c${chr_num}_b0_v1.bim" \
    -iin="${path_to_500kwes_plink_files}/ukb23158_c${chr_num}_b0_v1.fam" \
    -icmd="plink2 --bfile ukb23158_c${chr_num}_b0_v1 --score ${basename}.txt no-mean-imputation --no-fam-pheno --out ${basename}" \
  	--destination "/project/data/countfiles_perCHR_MAF0.01/" \
    --instance-type mem1_ssd1_v2_x16 --priority normal \
    --name "${basename}" --tag "count_rare_allCHR_${type}" \
    -y 
done
