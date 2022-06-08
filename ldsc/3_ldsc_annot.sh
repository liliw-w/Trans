#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --output=o.%x.%J.txt
#SBATCH --error=e.%x.%J.txt
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=liliw1@uchicago.edu
#SBATCH --mem=30G
#SBATCH --job-name=ldsc

#source activate ldsc
which python


prefix=$1

for chr in {1..22}
do
~/.conda/envs/ldsc/bin/python ldsc.py \
  --l2 \
  --bfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
  --ld-wind-cm 1 \
  --annot ldsc_annot/${prefix}.${chr}.annot.gz \
  --thin-annot \
  --out ldsc_annot/${prefix}.${chr} \
  --print-snps data/w_hm3.snplist.snponly
echo "Partitioned ldsc on chr"${chr}" is done."
done

# --print-snps only print LD Scores for the SNPs listed (one ID per row)
