#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --output=o.%x.%J.txt
#SBATCH --error=e.%x.%J.txt
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=liliw1@uchicago.edu
#SBATCH --mem=10G
#SBATCH --job-name=powerN

module load R/3.6.1
cd /home/liliw1/scratch/simulation_lambda0.1

Rscript --no-restore --no-save power_cal.R \
new_Sigma/simulation.alt.N.lambda0.1.varb1e-3.K101.rds \
new_Sigma/simulation.null.lambda0.1.K101.rds \
new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds \
0.1 \
True

