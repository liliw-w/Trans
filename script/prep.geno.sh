dir_plink='/project2/xuanyao/data/DGN/genotype/'
dir_data='/project2/xuanyao/llw/DGN/data/'
cd $dir_data

cp $dir_plink'chr'*'_QCed.'{bed,bim,fam} .

for chr in `seq 22`
do
cat 'chr'$chr'_QCed.bim' | \
awk '{key=sprintf("%s:%s", $1, $4); printf("chr%s\t%s\t%s\t%s\t%s\t%s\t\n", $1,key,$3,$4,$5,$6)}' > 'chr'$chr'_QCed_dup.bim'
rm 'chr'$chr'_QCed.bim'
mv 'chr'$chr'_QCed_dup.bim' 'chr'$chr'_QCed.bim'
echo chr$chr
done
