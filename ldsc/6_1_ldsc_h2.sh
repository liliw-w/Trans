#source activate ldsc
which python


prefix=$1

for gwasPhenocode in 30080 30090 30100 30110 30010 30020 30030 30040 30050 30060 30070 30270 30240 30250 30260 30280 30290 30300 30000 30120 30130 30140 30150 30160 30180 30190 30200 30210 30220
do
echo 'Trait '${gwasPhenocode}' Module'${prefix}'is running.'
~/.conda/envs/ldsc/bin/python ldsc.py \
  --h2 gwas/${gwasPhenocode}.sumstats.gz \
  --w-ld-chr data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --ref-ld-chr data/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,ldsc_annot/${prefix}. \
  --overlap-annot \
  --frqfile-chr data/1000G_Phase3_frq/1000G.EUR.QC. \
  --print-coefficients \
  --out h2_enrich_par/${gwasPhenocode}_${prefix}_baseline
done


#	--overlap-annot Define functional categories
#	--ref-ld-chr LD scores (baseline model LD scores)
#	--w-ld-chr regression weights
#	--frqfile-chr allele frequencies
#	--h2 GWAS summary statistics
#	--out output

