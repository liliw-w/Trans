#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --output=o.%x.%J.txt
#SBATCH --error=e.%x.%J.txt
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=liliw1@uchicago.edu
#SBATCH --mem=30G
#SBATCH --job-name=convgwas

#source activate ldsc
which python

trait=$1

~/.conda/envs/ldsc/bin/python munge_sumstats.py \
  --sumstats gwas/${trait}.tsv.gz \
  --merge-alleles data/w_hm3.snplist \
  --out gwas/${trait} \
  --a1-inc \
  --chunksize 500000

# --merge-alleles recommend only keeping HapMap3 SNPs

#1. add rs id to nealelab sum stats (by .bim file)
#2. add N
#3. check ref allele and second allele. If they are the same as in .bim and hapmap. (sign of z-scores)
#4. check how many hapmap snps (GRCh37) are in sum stats. If they align.


