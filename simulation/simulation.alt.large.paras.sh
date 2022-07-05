#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --output=o.%x.%J.txt
#SBATCH --error=e.%x.%J.txt
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=liliw1@uchicago.edu
#SBATCH --mem=30G
#SBATCH --job-name=large

module load R/3.6.1
cd /scratch/midway2/liliw1/tmp

Rscript --no-restore --no-save /home/liliw1/Trans/simulation/simulation.alt.large.paras.R
