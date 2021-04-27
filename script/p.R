#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#data.type = parameter[1]
#module.id = as.numeric(parameter[2])
#chr = as.numeric(parameter[3])
#p.method = parameter[4]

file.ex.var.regressed = snakemake@input[['file_ex_var_regressed']]
file.gene.meta = snakemake@input[['file_gene_meta']]
file.coexp.module = snakemake@input[['file_coexp_module']]
file.z = snakemake@input[['file_z']]
file.p = snakemake@output[['file_p']]
params1 = snakemake@params[['dir_script']]
chr = as.numeric(snakemake@params[['chr']])
module = as.numeric(snakemake@params[['module']])

library(data.table)

datExpr = readRDS(file.ex.var.regressed)
gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module = readRDS(file.coexp.module)

# zscore input
z.mat = fread(file.z)
z.mat = as.matrix(z.mat, rownames = TRUE)

# Apply new PCO
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))



gene_in_cluster = data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene.meta)
gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]


if(length(gene_trans) > 1){
  Sigma = cor(datExpr[, gene_trans])
  SigmaO = ModifiedSigmaOEstimate(Sigma) 

  z.mat_trans = z.mat[, gene_trans]; rm(z.mat)
  p.all = ModifiedPCOMerged(Z.mat=z.mat_trans,
                            Sigma=Sigma, SigmaO=SigmaO)
}else{
  p.all = NULL
}


saveRDS(p.all, file = file.p)

