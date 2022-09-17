############################################################
########### Check how many cross map genes are there for cis genes of each SNP ###########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


### I/O
file_snp_meta <- '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'
file_gene_anno <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
file_gene_meta <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/gene.meta.txt"
file_cross_map <- 'tmp_cross_map.txt.gz'
#'/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz'
source('~/Trans/followup/theme_my_pub.R')

file_num_cross_map_genes <- "num_cross_map_genes.rds"
file_cross_map_genes <- "cross_map_genes.rds"


snp_meta <- fread(file_snp_meta)
gene_anno <- fread(file_gene_anno)
gene_meta <- fread(file_gene_meta)
cross_map <- fread(file_cross_map, header = FALSE, col.names = c("gene1", "gene2", "score"))

#cross_map <- cross_map %>%
#  separate(gene1, c('gene1', NA), sep = "[.]") %>%
#  separate(gene2, c('gene2', NA), sep = "[.]")
gene_meta <- gene_meta %>% separate(GeneNameConv, c('GeneNameConv', NA), sep = "[.]")
gene_anno <- gene_anno %>% separate(Geneid, c('Geneid', NA), sep = "[.]")
cross_map <- cross_map %>% filter(gene1 %in% gene_meta$GeneNameConv | gene2 %in% gene_meta$GeneNameConv)
cross_map <- cross_map %>%
  left_join(gene_anno %>% select(Geneid, GeneSymbol),
                                     by = c("gene1" = "Geneid")) %>%
  rename("gene1_GeneSymbol" = "GeneSymbol") %>%
  left_join(gene_anno %>% select(Geneid, GeneSymbol),
            by = c("gene2" = "Geneid")) %>%
  rename("gene2_GeneSymbol" = "GeneSymbol")


### 1. cis- genes of each snp
dis_cis <- 1e+5
tmp <- NULL
for(i in 1:nrow(snp_meta)){
  x <- snp_meta[i, ]
  tmp_gene_anno <- gene_anno %>% filter(Chromosome == paste0("chr", x$SNPChr) &
                                          abs(Start - as.numeric(x$SNPPos))<dis_cis/2 )
  tmp <- rbind(tmp,
              cbind( paste(tmp_gene_anno$Geneid, collapse = ";"), paste(tmp_gene_anno$GeneSymbol, collapse = ";") )
              )
}
snp_meta_cis <- cbind(snp_meta, as_tibble(tmp) %>% setNames(c("cis_Geneid", "cis_GeneSymbol")) )
colnames(snp_meta_cis)
dim(snp_meta_cis)


### 2. cross map genes of cis genes of each SNP
res_num <- NULL
res_gene <- NULL
for(i in 1:nrow(snp_meta_cis)){
  x <- snp_meta_cis[i, ]
  cis_gene <- strsplit(x$cis_Geneid, split = ";")[[1]]
  tmp_cross_map <- cross_map %>% filter( (gene1 %in% cis_gene & gene2 %in% gene_meta$GeneNameConv) |
                                          (gene1 %in% gene_meta$GeneNameConv & gene2 %in% cis_gene) )
  
  cross_map_gene_uniq <- unique( c(tmp_cross_map$gene1, tmp_cross_map$gene2) )
  cross_map_gene_GeneSymbol_uniq <- unique( c(tmp_cross_map$gene1_GeneSymbol, tmp_cross_map$gene2_GeneSymbol) )
  
  
  str(cis_gene)
  str(cross_map_gene_uniq)
  print(sum(cross_map_gene_uniq %in% gene_meta$GeneNameConv))
  print(sum(cross_map_gene_uniq %in% cis_gene))
  res_num <- rbind(res_num,
              c(length(cis_gene), sum(cross_map_gene_uniq %in% gene_meta$GeneNameConv)))
  res_gene <- rbind(res_gene,
                    cbind(x,
                          "cross_Geneid" = paste(cross_map_gene_uniq, collapse = ";"),
                          "cross_GeneSymbol" = paste(cross_map_gene_GeneSymbol_uniq, collapse = ";")
                          )
                    )
  
  cat("SNP", i, "\n")
}
saveRDS(as_tibble(res_num), file_num_cross_map_genes)
saveRDS(res_gene, file_cross_map_genes)





