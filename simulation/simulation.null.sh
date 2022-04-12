#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --output=o.%x.%J.txt
#SBATCH --error=e.%x.%J.txt
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=liliw1@uchicago.edu
#SBATCH --mem=55G
#SBATCH --job-name=null

module load R/3.6.1
cd /scratch/midway2/liliw1/simulation_lambda0.1

Rscript --no-restore --no-save simulation.null.R
