signals_file=$1
sig_uniq=$2
sig_indp=$3
dir_geno=$4
geno_prefix=$5
geno_suffix=$6


module load plink
module load R/3.6.1

Rscript ./script/signals.uniq.R $signals_file $sig_uniq

for chr in {1..22}
do
plink --bfile $dir_geno$geno_prefix$chr$geno_suffix \
--extract $sig_uniq \
--indep-pairwise 50 5 0.2 \
--out indep.eigene.eqtls.$chr
echo chr$chr
done

cat indep.eigene.eqtls.*.prune.in > in.txt
cat indep.eigene.eqtls.*.prune.out > out.txt

grep -vFf out.txt $sig_uniq > $sig_indp

rm indep.eigene.eqtls* in.txt out.txt

echo $dir_geno
echo "(module, snp):"$(wc $signals_file -l)
echo "unique snps:"$(wc $sig_uniq -l)
echo "independent snps:"$(wc $sig_indp -l)
