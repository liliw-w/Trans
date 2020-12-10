rm(list = ls())

parameter = commandArgs(trailingOnly = T)
gz.exp = parameter[1]
out = parameter[2]

require(data.table)
data.gz.exp = fread(gz.exp, header = TRUE)
ex = as.matrix(data.gz.exp[, -(1:4)])
rownames(ex) = data.gz.exp$gene_id
ex = t(ex)

saveRDS(ex, file = out)
