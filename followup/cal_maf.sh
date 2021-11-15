# calculate MAF
# /project2/xuanyao/llw/DGN/data/

$dir_geno=$1
cd $dir_geno

for chr in {1..22}
do
plink --bfile chr${chr}_QCed --freq --out chr${chr}
done

head -n 1 chr1.frq > chr_all.frq && tail -n +2 -q chr{1..22}.frq >> chr_all.frq
