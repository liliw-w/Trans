#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --output=o.%x.%J.txt
#SBATCH --error=e.%x.%J.txt
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=liliw1@uchicago.edu
#SBATCH --mem=30G
#SBATCH --job-name=N

module load R/3.6.1
cd /scratch/midway2/liliw1/simulation_lambda0.1

Rscript --no-restore --no-save simulation.alt.R 'N' 0.1 './script_lambda0.1/' './new_Sigma/Sigma-new_DGN_module29_K101.rds' './new_Sigma/simulation.alt.N.lambda0.1.varb1e-3.K101.rds' 0.001 0.3 500
