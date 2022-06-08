#source activate ldsc
which python

prefix=$1

for chr in {1..22}
do
~/.conda/envs/ldsc/bin/python make_annot.py \
		--gene-set-file geneset/${prefix}.GeneSet \
		--gene-coord-file data/ENSG_coord.txt \
		--windowsize 500 \
		--bimfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
		--annot-file ldsc_annot/${prefix}.${chr}.annot.gz
echo "Annot on chr"${chr}" is done."
done

# --gene-set-file a gene set file with the names of the genes in your gene set, one line per gene name.
# --gene-coord-file a gene coordinate file, with columns GENE, CHR, START, and END
# --windowsize The annotation will include all SNPs within this many base pairs of the transcribed region.
# --bimfile the plink bim file of the dataset you will use to compute LD scores.
# --annot-file the name of the annot file to output

