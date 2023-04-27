#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#module = as.numeric(parameter[1])
#chr = as.numeric(parameter[2])
#p.method = parameter[3]
#fdr.level = 0.05/18/22

input1 = snakemake@input[['file_p']]
input2 = snakemake@input[['file_p_null']]
output1 = snakemake@output[['file_q']]
#fdr.level = as.numeric(snakemake@params[['fdr_level']])

str(input1)
str(input2)

require(data.table)

p.obs = rbindlist(lapply(input1, function(x)
{tmp_y=readRDS(x);
 tmp_y = tmp_y[tmp_y < 1e-5]
if(!is.null(tmp_y) & length(tmp_y) > 0){
  as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)
}
}))

p.null = rbindlist(lapply(input2, function(x)
{tmp_y=readRDS(x);
 tmp_y = tmp_y[tmp_y < 1e-5]
if(!is.null(tmp_y) & length(tmp_y) > 0){
  as.data.table(setNames(tmp_y, paste0("n", strsplit(x, '.', fixed = T)[[1]][3], ":", names(tmp_y))), keep.rownames=T)
}
}))


p.obs.rank = frank(p.obs, V2)
all.rank = frank(rbindlist(list(p.obs, p.null)), V2)
q = pmin(all.rank[1:length(p.obs.rank)]/p.obs.rank - 1, 1)

res = data.frame("snp" = p.obs$V1, "p" = p.obs$V2, "q" = q)

saveRDS(res, output1)

#fwrite(res[q<fdr.level, ], output1,
#       sep = "\t", row.names = FALSE, col.names = FALSE)
