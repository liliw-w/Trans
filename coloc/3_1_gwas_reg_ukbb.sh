#!/bin/sh
#SBATCH --time=5:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --output=gwas_coloc.o

Rscript --no-restore --no-save ~/Trans/coloc/3_1_gwas_reg_ukbb.R

