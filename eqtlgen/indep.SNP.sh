sig_uniq=$1
sig_indp=$2
dir_geno=/project2/xuanyao/llw/DGN_data/geno/
geno_prefix=chr
geno_suffix=_QCed

#module load plink
#module load R/3.6.1
#set -o nounset -o errexit -o pipefail

for chr in {1..22}
do
plink --bfile $dir_geno$geno_prefix$chr$geno_suffix \
--extract $sig_uniq \
--indep-pairwise 50 5 0.2 \
--out indep.eigene.eqtls.$chr
#echo chr$chr
done

cat indep.eigene.eqtls.*.prune.in > in.txt
cat indep.eigene.eqtls.*.prune.out > out.txt
grep -vFf out.txt $sig_uniq | sort > $sig_indp
rm indep.eigene.eqtls* in.txt out.txt

