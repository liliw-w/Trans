#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#module = as.numeric(parameter[1])
#chr = as.numeric(parameter[2])
#p.method = parameter[3]
#fdr.level = 0.05/18/22

input1 = snakemake@input[['file_p']]
input2 = snakemake@input[['file_p_null']]
output1 = snakemake@output[['file_q']]
output2 = snakemake@output[['file_signals']]
module = as.numeric(snakemake@params[['module']])
fdr.level = as.numeric(snakemake@params[['fdr_level']])/as.numeric(snakemake@params[['Nmodule']])/as.numeric(snakemake@params[['Nchr']])


p.obs = readRDS(input1)
p.null = readRDS(input2)

names(p.obs) = paste0("C", module, ":", names(p.obs))
p.obs.rank = rank(p.obs)
all.rank = rank(c(p.obs, p.null))
q = pmin(all.rank[names(p.obs)]/p.obs.rank[names(p.obs)]-1, 1)

res = data.frame("snp" = names(p.obs), "p" = p.obs, "q" = q[names(p.obs)])

write.table(res, output1,
            sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(res[res$q < fdr.level, ], output2,
            sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
