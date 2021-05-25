#!/bin/bash

#module load R/3.6.1
#module load python
#module load plink

set -o nounset -o errexit -o pipefail

~/.conda/envs/snakemake/bin/snakemake \
-s Snakefile0 \
-j
