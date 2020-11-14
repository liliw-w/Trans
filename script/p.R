#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#data.type = parameter[1]
#module.id = as.numeric(parameter[2])
#chr = as.numeric(parameter[3])
#p.method = parameter[4]

input1 = snakemake@input[['file_ex_var_regressed']]
input2 = snakemake@input[['file_gene_meta']]
input3 = snakemake@input[['file_coexp_module']]
input4 = snakemake@input[['file_z']]
output1 = snakemake@output[['file_p']]
params1 = snakemake@params[['dir_script']]
chr = as.numeric(snakemake@params[['chr']])
module = as.numeric(snakemake@params[['module']])

library(data.table)

datExpr = readRDS(input1)
gene.meta = read.table(input2, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F)
gene_module = readRDS(input3)

# zscore input
z.mat = fread(input4)
z.mat = as.matrix(z.mat, rownames = TRUE)

# Apply new PCO
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))



gene_in_cluster = data.frame("gene" = names(gene_module$moduleLabels[gene_module$moduleLabels == module]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene.meta, by.x = 1, by.y = 1)
gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]


Sigma = cor(datExpr[, gene_trans])
SigmaO = ModifiedSigmaOEstimate(Sigma) 

z.mat_trans = z.mat[, gene_trans]; rm(z.mat)
p.all = ModifiedPCOMerged(Z.mat=z.mat_trans,
                          Sigma=Sigma, SigmaO=SigmaO)


saveRDS(p.all, file = output1)

