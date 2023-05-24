module load plink


dir_geno=/project2/xuanyao/llw/DGN_data/geno/
geno_prefix=chr
geno_suffix=_QCed

fdr=10

file_uniq_signal=/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/LD.prun.in.chr.module.perm10.fdr${fdr}.txt
realpath /project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module*_fdr${fdr}.txt > msig_all_signal_module_fdr${fdr}.txt
tr '\n' , < msig_all_signal_module_fdr${fdr}.txt

clump_file='/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module10_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module11_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module12_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module13_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module14_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module15_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module16_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module17_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module18_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module19_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module21_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module22_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module23_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module24_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module25_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module27_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module28_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module29_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module2_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module31_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module32_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module34_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module35_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module37_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module38_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module39_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module3_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module40_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module41_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module42_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module43_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module44_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module45_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module46_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module47_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module49_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module4_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module5_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module7_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module8_fdr10.txt,/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/ld_clump/module9_fdr10.txt'


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

awk 'FNR>1' clumped_*.clumped | grep "\S" > msig_clumped_allchr_fdr${fdr}.txt
rm clumped_*.nosex
rm clumped_*.clumped
rm clumped_*.log

