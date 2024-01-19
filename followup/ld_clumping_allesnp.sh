module load plink


dir_geno=/project2/xuanyao/llw/DGN_data/geno/
geno_prefix=chr
geno_suffix=_QCed


file_uniq_signal=/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/ld_clump/uniq_snp.txt

clump_file='/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/ld_clump/ld_clump_esnp.txt'

file_clump_res=clumped_esnp_allchr.txt


for chr in {1..22}
do
plink --bfile $dir_geno$geno_prefix$chr$geno_suffix \
      --extract $file_uniq_signal \
         --clump ${clump_file} \
         --clump-p1 1 \
         --clump-p2 1 \
         --clump-r2 0.2 \
         --clump-kb 250  \
        --out clumped_${chr}
#echo chr$chr
done

awk 'FNR>1' clumped_*.clumped | grep "\S" > $file_clump_res
rm -rf clumped_*.nosex
rm -rf clumped_*.clumped
rm -rf clumped_*.log

