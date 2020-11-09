rm(list = ls())

parameter = commandArgs(trailingOnly = T)
data.type = parameter[1]
module.id = as.numeric(parameter[2])
chr = as.numeric(parameter[3])
p.method = parameter[4]


if(data.type == "obs"){
  data.type = ""
}else if (data.type == "null"){
  data.type = paste0(data.type, ".")
}else{stop("Please specify a valid data.type: obs, null.\n")}


datExpr = readRDS("./data/ex_var_regressed.rds")
gene.meta = read.table("./data/gene.meta.txt", sep = '\t', header = T, stringsAsFactors = F)
gene_module = readRDS("./data/coexp.module.rds")

# zscore input
z.mat = fread(paste0("./z/z.", data.type, "module", module_id, ".chr", chr, ".txt.gz"))
z.mat = as.matrix(z.mat, rownames = TRUE)

# Apply new PCO
source(paste0("./script/ModifiedPCOMerged.R"))
source(paste0("./script/liu.r"))
source(paste0("./script/liumod.R"))
source(paste0("./script/davies.R"))
dyn.load("./script/qfc.so")
source(paste0("./script/ModifiedSigmaOEstimate.R"))



gene_in_cluster = data.frame("gene_name" = names(gene_module$moduleLabels[gene_module$moduleLabels == module.id]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene.meta, by.x = 1, by.y = 1)
gene_trans = gene_w_pos[gene_w_pos[,2] != paste0("chr", chr), "gene_name"]


Sigma = cor(datExpr[, gene_trans])
SigmaO = ModifiedSigmaOEstimate(Sigma, p.method = p.method) 

z.mat_trans = z.mat[, gene_trans]; rm(z.mat)
p.all = ModifiedPCOMerged(Z.mat=z.mat_trans,
                          Sigma=Sigma, SigmaO=SigmaO, p.method = p.method)


saveRDS(p.all, file = paste0("./p/p.", data.type, "module", module.id, ".chr", chr, ".", p.method, ".rds"))

