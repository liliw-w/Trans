rm(list = ls())

parameter = commandArgs(trailingOnly = T)
data.type = parameter[1]
module.id = parameter[2]
chr = parameter[3]

if(data.type == "obs"){
  data.type = ""
}else if (data.type == "null"){
  data.type = paste0(data.type, ".")
}else{stop("Please specify a valid data.type: obs, null.\n")}


library(data.table)
z = fread(paste0("./z/module", module.id, ".chr", chr, ".trans_qtl_pairs_z.txt"),
          header = TRUE, sep = "\t")
z.mat = dcast(z, snp~gene, value.var = "zscore", fun.aggregate = max, drop = FALSE)
fwrite(z.mat, paste0("./z./z.", data.type, "module", module.id, ".chr", chr, ".txt.gz"),
       sep = "\t", row.names = FALSE, col.names = TRUE)
