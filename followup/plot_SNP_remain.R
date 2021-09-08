rm(list = ls())
library(data.table)
library(tidyverse)
require(ggplot2)

file_gene_annotation = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
file_cross_mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"
file.ex.var.regressed = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/ex.var.regressed.rds"
file.gene.meta = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/gene.meta.txt"
file.coexp.module = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/coexp.module.rds"
#file.z = "/scratch/midway2/liliw1/DGN_PCO.lambda.01_real/z/z.module68.chr22.txt.gz"
file.p = NULL
params1 = "/scratch/midway2/liliw1/DGN_PCO.lambda.01_real/script/"


datExpr = readRDS(file.ex.var.regressed)
gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module = readRDS(file.coexp.module)
gene_annotation = fread(file_gene_annotation, header = TRUE)
cross_map = fread(file_cross_mappability)

seq.module = c(1, 10, 30, 50, 100, 130, 157, 158)
seq.chr = 22:1
mat_SNP_chr = matrix(nrow = length(seq.module), ncol = length(seq.chr),
                       dimnames = list(paste0("module", seq.module), paste0("chr", seq.chr)))
mat_SNP_remain = matrix(nrow = length(seq.module), ncol = length(seq.chr),
                 dimnames = list(paste0("module", seq.module), paste0("chr", seq.chr)))

for(chr in seq.chr){
  file.z = paste0("/scratch/midway2/liliw1/DGN_PCO.lambda.01_real/z/z.module68.chr", chr, ".txt.gz")
  z.mat = fread(file.z)
  z.mat = as.matrix(z.mat, rownames = TRUE)

  df_SNP <- data.frame("SNP" = rownames(z.mat),
                       "pos" = sapply(rownames(z.mat), function(x) as.numeric(strsplit(x, ":", fixed = FALSE)[[1]][2])),
                       check.names = FALSE, stringsAsFactors = FALSE)
  gene_annotation_chr = gene_annotation[gene_annotation$Chromosome == paste0("chr", chr), ]
  cis_Genes = lapply(df_SNP$pos, function(x) gene_annotation_chr[abs(x - gene_annotation_chr$Start)<(5*10^5), Geneid] )

  for(module in seq.module){
    gene_in_cluster = data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
    gene_w_pos = merge(gene_in_cluster, gene.meta)
    gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]


    module_Geneid = gene_annotation[match(gene_trans, gene_annotation$GeneSymbol), Geneid]
    ind_cross_map = (cross_map$V1 %in% module_Geneid | cross_map$V2 %in% module_Geneid) & (cross_map$V1 %in% gene_annotation_chr$Geneid | cross_map$V2 %in% gene_annotation_chr$Geneid)
    cross_map_compare = cross_map[ind_cross_map, ]
    cross_map_gene_chr = unique(c(cross_map_compare[cross_map_compare$V1 %in% gene_annotation_chr$Geneid, V1], cross_map_compare[cross_map_compare$V2 %in% gene_annotation_chr$Geneid, V2]))

    SNP_remain <- sapply(cis_Genes, function(x){
      sum(x %in% cross_map_gene_chr) == 0
    } )

    mat_SNP_remain[paste0("module", module), paste0("chr", chr)] = sum(SNP_remain)
    mat_SNP_chr[paste0("module", module), paste0("chr", chr)] = length(SNP_remain)

    cat("Module ", module, " is done!", "\n")
  }
  cat("Chr ", chr, " is done!", "\n")
}
save(mat_SNP_chr, mat_SNP_remain, file = "SNP_remain.RData")


sapply(cis_Genes, function(x) length(x) )
SNP_remain2 <- sapply(cis_Genes, function(x){
  sum(x %in% cross_map_gene_chr)
} )

png("a.png")
hist(sapply(cis_Genes, function(x) length(x) ))
dev.off()

png("b.png")
hist(SNP_remain2)
dev.off()

### Prepare plot data
df_SNP_remove = cbind("Module" = rownames(mat_SNP_chr), as.data.frame(mat_SNP_chr-mat_SNP_remain, check.names = FALSE, stringsAsFactors = FALSE)) %>%
  pivot_longer(!Module, names_to = "Chr", values_to = "SNP")
df_SNP_remove$type = "SNP_remove"

df_SNP_remain = cbind("Module" = rownames(mat_SNP_remain), as.data.frame(mat_SNP_remain, check.names = FALSE, stringsAsFactors = FALSE)) %>%
  pivot_longer(!Module, names_to = "Chr", values_to = "SNP")
df_SNP_remain$type = "SNP_remain"

df_SNP_chr = cbind("Module" = rownames(mat_SNP_chr), as.data.frame(mat_SNP_chr, check.names = FALSE, stringsAsFactors = FALSE)) %>%
  pivot_longer(!Module, names_to = "Chr", values_to = "SNP")
df_SNP_chr$type = "SNP_chr"


df_SNP_num = rbind(df_SNP_remain, df_SNP_remove)
df_SNP_num$Chr = factor(df_SNP_num$Chr, levels = paste0("chr", 1:22))
df_SNP_num$Module = factor(df_SNP_num$Module, levels = paste0("module", seq.module))

df_SNP_percent = merge(df_SNP_remain, df_SNP_chr, by = c("Module", "Chr"), suffixes = c("_remain", "_chr"))
df_SNP_percent$Chr = factor(df_SNP_percent$Chr, levels = paste0("chr", 1:22))
df_SNP_percent$Module = factor(df_SNP_percent$Module, levels = paste0("module", seq.module))


### Plot: #SNPs remained & removed on each Chr for each module
fig_bar <- ggplot(df_SNP_num, aes(Chr, SNP)) +
  facet_wrap(~Module) +
  geom_bar(aes(fill = type), stat="identity") +
  labs(y = "#SNPs", fill = NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top") +
  scale_fill_manual(values=c("dodgerblue2", "tomato3"))
ggsave("plot/SNP_remain_bar.png", fig_bar)


### Plot: Percentage of SNPs remained on each Chr for each module
fig_percent <- ggplot(df_SNP_percent, aes(Chr, SNP_remain/SNP_chr, group = Module, color = Module)) +
  geom_line() +
  labs(y = "% SNPs", color = NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top")
ggsave("plot/SNP_remain_percent.png", fig_percent)


system("bash ~/imgcat plot/SNP_remain_bar.png plot/SNP_remain_percent.png")
