##############################################
########### ld clumping for gwas ###########
##############################################
set -e # exit on first error

module load plink

# gwas file list
pheno_code_seq=(30080 30090 30100 30110 30010 30020 30030 30040 30050 30060 30070 30270 30240 30250 30260 30280 30290 30300 30000 30120 30130 30140 30150 30160 30180 30190 30200 30210 30220)


# gwas file
file_uniq_signal_prefix='ld_clump/uniq_snp_ukb_nealelab'
clump_file_prefix='ld_clump/ld_clump_snp_withp'


# geno file
dir_geno='/project2/xuanyao/llw/DGN_data/geno/'
geno_prefix='chr'
geno_suffix='_QCed_filtdup'


for gwasPhenocode in "${pheno_code_seq[@]}"
do
  # trait info
  gwas_trait_type='continuous'
  echo $gwasPhenocode
  
  # output file
  file_clump_res=ld_clump/clumped_snp_allchr_${gwas_trait_type}-${gwasPhenocode}.txt
  
  
  # ld clumping
  for chr in {1..22}
  do
  plink --bfile $dir_geno$geno_prefix$chr$geno_suffix \
        --extract ${file_uniq_signal_prefix}_chr${chr}.txt \
        --clump ${clump_file_prefix}_${gwas_trait_type}-${gwasPhenocode}.txt.gz \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.2 \
        --clump-kb 250  \
        --out clumped_${gwas_trait_type}-${gwasPhenocode}_${chr}
  #echo chr$chr
  done
  
  # merge file
  awk 'FNR>1' clumped_${gwas_trait_type}-${gwasPhenocode}_*.clumped | grep "\S" > $file_clump_res
  rm -rf clumped_${gwas_trait_type}-${gwasPhenocode}_*.nosex
  rm -rf clumped_${gwas_trait_type}-${gwasPhenocode}_*.clumped
  rm -rf clumped_${gwas_trait_type}-${gwasPhenocode}_*.log
done

