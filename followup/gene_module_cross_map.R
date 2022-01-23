rm(list = ls())
library(data.table)
library(tidyverse)

file_cross_map = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"
file_map = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_gene_mappability.txt.gz"
file_gene_meta = "/scratch/midway2/liliw1/DGN_no_filter_on_mappability/result/gene.meta.txt"
file_coexp = "/scratch/midway2/liliw1/DGN_no_filter_on_mappability/result/coexp.module.rds"

module = 51
gene_of_interest = "SENP7"

cross_map = fread(file_cross_map, header = FALSE, col.names = c("gene1", "gene2", "score"))
map = fread(file_map, header = FALSE, col.names = c("gene", "score_map"))
gene_meta = fread(file_gene_meta, header = TRUE)
coexp = readRDS(file_coexp)$moduleLabels


gene_in_module = names(coexp[coexp == module])
gene_meta_in_module = gene_meta %>% filter(gene %in% gene_in_module)
gene_meta_of_interest = gene_meta %>% filter(gene %in% gene_of_interest)

tmp = unique(c(cross_map$gene1, cross_map$gene2))
sum(gene_meta_in_module$GeneNameConv %in% tmp)
sum(gene_meta_of_interest$GeneNameConv %in% tmp)

cross_map_of_interest = cross_map %>%
  filter(gene1 %in% gene_meta_of_interest$GeneNameConv | gene2 %in% gene_meta_of_interest$GeneNameConv) %>%
  filter(gene1 %in% gene_meta_in_module$GeneNameConv | gene2 %in% gene_meta_in_module$GeneNameConv)

ind_filp = cross_map_of_interest$gene1 != gene_meta_of_interest$GeneNameConv
cross_map_of_interest[ind_filp, 2] = cross_map_of_interest[ind_filp, 1]
cross_map_of_interest$gene1 = gene_meta_of_interest$GeneNameConv

cross_map_of_interest = cross_map_of_interest %>%
  left_join(gene_meta %>% select(c(gene, chr, GeneNameConv)), by = c("gene1" = "GeneNameConv")) %>%
  rename("gene1_name" = "gene", "chr1" = "chr") %>%
  left_join(gene_meta %>% select(c(gene, chr, GeneNameConv)), by = c("gene2" = "GeneNameConv"), suffix = c("", "2") ) %>%
  rename("gene2_name" = "gene", "chr2" = "chr") %>%
  arrange(desc(score))

cross_map_of_interest = cross_map_of_interest %>% left_join(map, by = c("gene1" = "gene")) %>%
  rename("score_map1" = "score_map") %>%
  left_join(map, by = c("gene2" = "gene")) %>%
  rename("score_map2" = "score_map")
fwrite(cross_map_of_interest,
       file = paste0(gene_of_interest, "_module", module, "_size", length(gene_in_module), "_crossMap.txt"),
       quote = FALSE, sep = "\t")

cross_map_of_interest

fwrite(data.table(gene_in_module), "tmp.txt")

