#!/bin/bash
#SBATCH --job-name=coloc.29traits
#SBATCH --partition=broadwl
#SBATCH --output=coloc.%J
#SBATCH --error=coloc.%J
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=50G

#module load R/3.6.1
cd /scratch/midway2/liliw1/coloc/

Rscript --no-restore --no-save script/coloc.29traits.R
