#!/bin/bash
#SBATCH --job-name=qtl_coloc
#SBATCH --partition=bigmem2
#SBATCH --output=%x.%J
#SBATCH --error=%x.%J
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=200G

#module load R/3.6.1
cd /scratch/midway2/liliw1/coloc_MSigDB

Rscript --no-restore --no-save /home/liliw1/Trans/coloc/1_qtl_reg.R
