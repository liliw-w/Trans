rm(list = ls())

parameter = commandArgs(trailingOnly = T)
inp = parameter[1]
outp = parameter[2]

library(data.table)
z = fread(inp,
          header = TRUE, sep = "\t")
z.mat = dcast(z, snp~gene, value.var = "zscore", fun.aggregate = max, drop = FALSE)
fwrite(z.mat, outp,
       sep = "\t", row.names = FALSE, col.names = TRUE)
