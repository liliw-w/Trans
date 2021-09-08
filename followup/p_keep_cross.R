rm(list = ls())
library(data.table)

file_gene_annotation = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
file_cross_mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"
file.ex.var.regressed = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/ex.var.regressed.rds"
file.gene.meta = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/gene.meta.txt"
file.coexp.module = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/coexp.module.rds"
#file.z = "/scratch/midway2/liliw1/DGN_PCO.lambda.01_real/z/z.module68.chr22.txt.gz"
file.p = NULL
params1 = "/scratch/midway2/liliw1/DGN_PCO.lambda.01_real/script/"
chr = 22
module = 68


datExpr = readRDS(file.ex.var.regressed)
gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module = readRDS(file.coexp.module)
gene_annotation = fread(file_gene_annotation, header = TRUE)
cross_map = fread(file_cross_mappability)

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
gene_annotation_chr = gene_annotation[gene_annotation$Chromosome == paste0("chr", chr), ]


module_Geneid = gene_annotation[match(gene_trans, gene_annotation$GeneSymbol), Geneid]
ind_cross_map = (cross_map$V1 %in% module_Geneid | cross_map$V2 %in% module_Geneid) & (cross_map$V1 %in% gene_annotation_chr$Geneid | cross_map$V2 %in% gene_annotation_chr$Geneid)
cross_map = cross_map[ind_cross_map, ]
cross_map_gene_chr = unique(c(cross_map[cross_map$V1 %in% gene_annotation_chr$Geneid, V1], cross_map[cross_map$V2 %in% gene_annotation_chr$Geneid, V2]))



df_SNP <- data.frame("SNP" = rownames(z.mat),
                     "pos" = sapply(rownames(z.mat), function(x) as.numeric(strsplit(x, ":", fixed = FALSE)[[1]][2])),
                     check.names = FALSE, stringsAsFactors = FALSE)
cis_Genes = lapply(df_SNP$pos, function(x) gene_annotation_chr[abs(x - gene_annotation_chr$Start)<(5*10^5), Geneid] )


system.time(
  tmp <- sapply(cis_Genes[1:100], function(x){
    colSums(matrix(outer(x, module_Geneid, FUN = "paste", sep=";") %in% cross_map$V4 |
                     t(outer(module_Geneid, x, FUN = "paste", sep=";")) %in% cross_map$V4, ncol = length(module_Geneid)) ) == 0

  } )

)
tmp = colSums(tmp)
sum(tmp>0)

system.time(
  tmp2 <- sapply(cis_Genes[1:100], function(x){
    sum(x %in% cross_map_gene_chr)
  } )

)
sum(tmp2>0)


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

