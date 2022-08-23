#source activate ldsc
which python


gwasPhenocode='28067908_cd'

for prefix in M{1..166}
do
echo 'Trait '${gwasPhenocode}' Module'${prefix}' is running.'
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

