~/.conda/envs/snakemake/bin/snakemake \
--cluster-config /project2/xuanyao/llw/DGN_PCO.lambda.01/cluster_config.json \
--cluster "sbatch -J {cluster.job-name} -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} --mem={cluster.mem} -o {cluster.output} -e {cluster.error} --mail-type={cluster.email-type} --mail-user={cluster.email}" \
-s Snakefile_DGN_est_Sigma \
--jobs 95
#--dryrun
