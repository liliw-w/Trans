##############################################
########### remove dup snps from bfiles ###########
##############################################
set -e # exit on first error

module load plink


# geno file
dir_geno='/project2/xuanyao/llw/DGN_data/geno/'
geno_prefix='chr'
geno_suffix='_QCed_filtdup'


# remove dup snps
for chr in {1..22}
do
  cut -f 2 $dir_geno$geno_prefix$chr$geno_suffix.bim | sort | uniq -d > $geno_prefix$chr$geno_suffix.dups
  
  plink --bfile $dir_geno$geno_prefix$chr$geno_suffix \
        --exclude $geno_prefix$chr$geno_suffix.dups \
        --make-bed \
        --out ${dir_geno}${geno_prefix}${chr}${geno_suffix}_filtdup
  
  rm $geno_prefix$chr$geno_suffix.dups
done

