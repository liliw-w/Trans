rm(list = ls())
library(data.table)
library(tidyverse)
require(ggplot2)

### File I/O
file_gene_annotation = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
file_cross_mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"
file.ex.var.regressed = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/ex.var.regressed.rds"
file_out = "crossGene_gene.rds"

datExpr = readRDS(file.ex.var.regressed)
gene_annotation = fread(file_gene_annotation, header = TRUE)
cross_map = fread(file_cross_mappability)

gene_all = gene_annotation[match(colnames(datExpr), gene_annotation$GeneSymbol), Geneid]

crossGene_gene = rep(NA, length(gene_all)); names(crossGene_gene) = gene_all
for(k in 1:length(gene_all)){
  x = gene_all[k]
  crossGene_gene[k] = sum(cross_map$V1 %in% x | cross_map$V2 %in% x)

  cat(k, "-th gene searched!", "\n")
}
saveRDS(crossGene_gene, file_out)


### Plot: #cross genes for each individual genes across all genes in DGN
#crossGene_gene = readRDS("~/Downloads/crossGene_gene.rds")

#### Prepare data
df_crossGene_gene = cbind("Gene" = names(crossGene_gene), as.data.frame(crossGene_gene, check.names = FALSE, stringsAsFactors = FALSE))

group_breaks = sort(c(0, 1, 10, 100, 500, 1000, 2000, 3000, 4000, 5000, max(df_crossGene_gene$crossGene_gene)))
x_lable = paste0(c("=", "=", rep("<=", length(group_breaks)-2)), group_breaks)
df_crossGene_gene$group = factor(sapply(df_crossGene_gene$crossGene_gene, function(x)
  group_breaks[x <=  group_breaks][1] ), levels = group_breaks)

df_summary = df_crossGene_gene %>% group_by(group) %>% count()
df_summary$percent = round(df_summary$n / sum(df_summary$n)*100, digits = 0)

#### Plot
fig_crossGene_gene <-ggplot(df_crossGene_gene, aes(x = group)) +
  geom_bar() +
  geom_text(data = df_summary, aes(x = group, y = n+150, label = paste0(n,"\n","(",percent,"%)")) ) +
  scale_x_discrete(labels = x_lable) +
  labs(x = "#cross genes", y = "#genes") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top")

ggsave("plot/crossGene_gene.png", fig_crossGene_gene)
system("bash ~/imgcat plot/crossGene_gene.png")

