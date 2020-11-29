input1 = snakemake@input[['file_q']]
output2 = snakemake@output[['file_signals']]
fdr.level = as.numeric(snakemake@params[['fdr_thre']])

require(dplyr)

p.obs = lapply(input1, function(x) readRDS(x) )
res = bind_rows(p.obs) %>% group_by(snp, p) %>% summarise("q" = mean(q)) %>% filter(q < fdr.level)
write.table(res, output2,
            sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
