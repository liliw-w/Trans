#!/bin/bash

conda activate /scratch/midway2/liliw1/conda_env/rstudio-server

file_cluster_config=/home/liliw1/Trans/modules/cluster_config_inflation.json
file_snmk=/home/liliw1/Trans/modules/Snakefile_inflation


~/.conda/envs/snakemake/bin/snakemake \
--cluster-config ${file_cluster_config} \
--cluster "sbatch -J {cluster.job-name} -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} --mem={cluster.mem} -o {cluster.output} -e {cluster.error} --mail-type={cluster.email-type} --mail-user={cluster.email}" \
--jobs 95 \
-s ${file_snmk}
