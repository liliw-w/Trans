cd /scratch/midway2/liliw1/UKB_nealelab/

for i in 30060 30070 30270 30240 30250 30260 30280 30290 30300 30000 30120 30130 30140 30150 30160 30180 30190 30200 30210 30220
do
echo $i
zcat phenotype_manifest.tsv.bgz | cut -f 1-3,5,43 | grep $i | wc
cmd=$(zcat phenotype_manifest.tsv.bgz | cut -f 43 | grep $i)
eval ${cmd}
done
