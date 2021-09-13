rm(list = ls())
library(data.table)
library(tidyverse)
require(ggplot2)

file_gene_annotation = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
file_cross_mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"
file.gene.meta = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/gene.meta.txt"
file.coexp.module = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/coexp.module.rds"
file_out = "crossGene.RData"

seq.module = c(1, 5, 12, 40, 111, 112)

gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module = readRDS(file.coexp.module)
gene_annotation = fread(file_gene_annotation, header = TRUE)
cross_map = fread(file_cross_mappability)


seq.chr = 22:1
mat_crossGene_module = matrix(nrow = length(seq.module), ncol = 2,
                              dimnames = list(paste0("module", seq.module), c("crossGene", "moduleSize") ))
mat_crossGene_chr = matrix(nrow = length(seq.module), ncol = length(seq.chr),
                           dimnames = list(paste0("module", seq.module), paste0("chr", seq.chr)))
for(module in seq.module){
  gene_in_cluster = data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
  gene_w_pos = merge(gene_in_cluster, gene.meta)

  module_Geneid = gene_annotation[match(gene_w_pos[, "gene"], gene_annotation$GeneSymbol), Geneid]
  mat_crossGene_module[paste0("module", module), c("crossGene", "moduleSize")] =
    c(length(unique(c(cross_map[cross_map$V1 %in% module_Geneid, V2], cross_map[cross_map$V2 %in% module_Geneid, V1]))),
      length(module_Geneid))

  for(chr in seq.chr){
    gene_annotation_chr = gene_annotation[gene_annotation$Chromosome == paste0("chr", chr), ]

    gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]
    module_Geneid = gene_annotation[match(gene_trans, gene_annotation$GeneSymbol), Geneid]
    ind_cross_map = (cross_map$V1 %in% module_Geneid | cross_map$V2 %in% module_Geneid) & (cross_map$V1 %in% gene_annotation_chr$Geneid | cross_map$V2 %in% gene_annotation_chr$Geneid)
    cross_map_compare = cross_map[ind_cross_map, ]
    cross_map_gene_chr = unique(c(cross_map_compare[cross_map_compare$V1 %in% gene_annotation_chr$Geneid, V1], cross_map_compare[cross_map_compare$V2 %in% gene_annotation_chr$Geneid, V2]))

    ind_crossmodule = cross_map$V1 %in% module_Geneid | cross_map$V2 %in% module_Geneid
    mat_crossGene_chr[paste0("module", module), paste0("chr", chr)] =
      length(unique(c(cross_map[ind_crossmodule & (cross_map$V1 %in% gene_annotation_chr$Geneid), V1], cross_map[ind_crossmodule & (cross_map$V2 %in% gene_annotation_chr$Geneid), V2] )))

    cat("Chr ", chr, " is done!", "\n")
  }

  cat("Module ", module, " is done!", "\n")
  save(mat_crossGene_module, mat_crossGene_chr, file = file_out)
}


### Plot: #cross genes for each module
#### Prepare data
df_crossGene_module = cbind("Module" = rownames(mat_crossGene_module), as.data.frame(mat_crossGene_module, check.names = FALSE, stringsAsFactors = FALSE))
df_crossGene_module$Module = factor(df_crossGene_module$Module, levels = paste0("module", sort(seq.module)))

#### Plot
fig_crossGene_module <- ggplot(df_crossGene_module, aes(x = Module)) +
  geom_bar(aes(y = -moduleSize), stat="identity", fill = "grey65", position = position_nudge(y = -100)) +
  geom_text(aes(y = -moduleSize-100, label = moduleSize)) +
  geom_bar(aes(y = crossGene), stat="identity", fill = "grey51") +
  labs(x = "Module", y = "Cross genes") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top")
ggsave("plot/crossGene_module.png", fig_crossGene_module)
system("bash ~/imgcat plot/crossGene_module.png")


### Plot: #cross genes for each module on each chr
#### Prepare data
df_crossGene_chr = cbind("Module" = rownames(mat_crossGene_chr), as.data.frame(mat_crossGene_chr, check.names = FALSE, stringsAsFactors = FALSE)) %>%
  pivot_longer(!Module, names_to = "Chr", values_to = "crossGene_chr")
df_crossGene_chr$Chr = factor(df_crossGene_chr$Chr, levels = paste0("chr", 1:22))
df_crossGene_chr$Module = factor(df_crossGene_chr$Module, levels = paste0("module", sort(seq.module)))

#### Plot
fig_crossGene_chr <- ggplot(df_crossGene_chr, aes(Chr, crossGene_chr, group = Module, color = Module)) +
  geom_line() +
  labs(x = "Chr", y = "Cross genes on chr", color = NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top")
ggsave("plot/crossGene_chr.png", fig_crossGene_chr)
system("bash ~/imgcat plot/crossGene_chr.png")
