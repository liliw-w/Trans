#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#module = as.numeric(parameter[1])
#chr = as.numeric(parameter[2])
#p.method = parameter[3]
#fdr.level = 0.05/18/22

input1 = snakemake@input[['file_p']]
input2 = snakemake@input[['file_p_null']]
output1 = snakemake@output[['file_signals']]
module = as.numeric(snakemake@params[['module']])
fdr.level = as.numeric(snakemake@params[['fdr_level']])/as.numeric(snakemake@params[['Nmodule']])

str(input1)
str(input2)

require(data.table)

p.obs = rbindlist(lapply(input1, function(x) 
  {tmp_y=readRDS(x);
  as.data.table(setNames(tmp_y, paste0("C", module, ":", names(tmp_y))), keep.rownames=T)}))
p.null = rbindlist(lapply(input2, function(x) 
  {tmp_y=readRDS(x);as.data.table(tmp_y, keep.rownames=T)}))

p.obs.rank = frank(p.obs, V2)
names(p.obs.rank) = p.obs$V1

all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
names(all.rank) = c(p.obs$V1, p.null$V1)

q = pmin(all.rank[names(p.obs.rank)]/p.obs.rank-1, 1)

res = data.table("snp" = p.obs$V1, "p" = p.obs$V2, "q" = q[p.obs$V1])

fwrite(res[q<fdr.level, ], output1,
       sep = "\t", row.names = FALSE, col.names = FALSE)
