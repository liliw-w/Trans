#!/bin/bash

#module load R/3.6.1
#module load python
#module load plink

#set -o nounset -o errexit -o pipefail

dir_env_snkm=/scratch/midway3/liliw1/conda_env/snakemake


if [[ ! -d "logs/" ]]
    then
      mkdir logs
fi

${dir_env_snkm}/bin/snakemake \
-s ~/Trans/simulation/Snakefile \
--cluster-config ~/Trans/simulation/cluster_config.json \
--cluster "sbatch -p {cluster.partition} --account {cluster.account} --qos {cluster.qos} -J {cluster.job-name} -t {cluster.time} -N {cluster.nodes} --cpus-per-task={cluster.cpus-per-task} --mem={cluster.mem} -o {cluster.output} -e {cluster.error} --mail-type={cluster.email-type} --mail-user={cluster.email}" \
--jobs 95
#--allowed-rules add_cpeak_info aggregate_top_p plt_peak_group_distribution plt_top_p \

