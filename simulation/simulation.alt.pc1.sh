#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --output=o.%x.%J.txt
#SBATCH --error=e.%x.%J.txt
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=liliw1@uchicago.edu
#SBATCH --mem=30G
#SBATCH --job-name=pc1

module load R/3.6.1
cd /scratch/midway2/liliw1/tmp

Rscript --no-restore --no-save /home/liliw1/Trans/simulation/simulation.alt.pc1.R
