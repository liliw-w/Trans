rm(list = ls())

#parameter = commandArgs(trailingOnly = T)
#dir.data = parameter[1]
#data.type = parameter[2]

data.type = snakemake@params[['data_type']]
if(data.type == "obs"){
  data.type = ""
}else if (data.type == "null"){
  data.type = paste0(data.type, ".")
}else{stop("Please specify a valid data.type: obs, null.\n")}

# sample name
sample.name = read.table(paste0(dir.data, "sample.name.txt"), header=F, stringsAsFactors=F)

# expression
load(paste0(dir.data, "ex.rdata"))
if(data.type == "null."){
  perm = sample(nrow(ex))
  ex = ex[perm, ]
  
  cov_all = read.table(paste0(dir.data, "covariates.txt"),
                       sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cov_all.perm = cov_all[, perm]; colnames(cov_all.perm) = colnames(cov_all)
  write.table(cov_all.perm, paste0(dir.data, 'covariates.null.txt'), quote = F, sep = "\t")
}
colnames(ex) = gnames; rownames(ex) = sample.name[, 1]

# gene position
gene.pos = as.matrix(read.table(paste0(dir.data, "gene.meta.txt"), header=T, stringsAsFactors=F))
rownames(gene.pos) = gene.pos[, 1]
gene.pos = gene.pos[, c("chr", "start", "end", "gene")]
gene.pos[, "chr"] = paste0(gene.pos[, "chr"], "NA")
colnames(gene.pos)[1] = "#chr"

# module
module.info = readRDS("./data/coexp.module.rds")
Nmodule = max(module.info$moduleLabels)
for(k in 1:Nmodule){
  gene.in.module = names(module.info$moduleLabels)[module.info$moduleLabels==k]
  res = cbind(gene.pos[gene.in.module, ], t(ex[, gene.in.module]))
  
  write.table(res, paste0(dir.data, "expression.", data.type, "module",k,".bed"),
              quote = F, sep = "\t", row.name = F)
}
