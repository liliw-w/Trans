module load R/3.6.1

set -o nounset -o errexit -o pipefail

~/.conda/envs/snakemake/bin/snakemake \
--cluster-config cluster_config.json \
--cluster "sbatch -J {cluster.job-name} -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} --mem={cluster.mem} -o {cluster.output} -e {cluster.error} --mail-type={cluster.email-type} --mail-user={cluster.email}" \
--jobs 90
