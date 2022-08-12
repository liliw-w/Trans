rm(list = ls())
library(data.table)
library(tidyverse)

file.gene.meta = '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
file_signal_cis_genes = '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/signal_cis_genes.txt'
file_coexp_module <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds'
dis_cis <- 1e+6

gene.meta = fread(file.gene.meta, header = TRUE)
signal_cis_genes <- fread(file_signal_cis_genes, header = TRUE)
coexp_module <- readRDS(file_coexp_module)$moduleLabels

coexp_module <- enframe(coexp_module, "gene", "m")
gene.meta[gene.meta$Strand == "-", ]$Start <- gene.meta[gene.meta$Strand == "-", ]$End

cis_gene_module <- sapply(
  1:nrow(signal_cis_genes),
  function(x){
    tmp_cis = signal_cis_genes[x, ]
    tmp_coexp_module = coexp_module %>%
      filter(paste0("module", m) == tmp_cis$module) %>%
      left_join(
        gene.meta %>% select(GeneSymbol:Class) %>% distinct(GeneSymbol, .keep_all = TRUE),
        by = c("gene" = "GeneSymbol")
      ) %>%
      filter(Chromosome == paste0("chr", tmp_cis$chr) ) %>%
      mutate("dis" = abs(tmp_cis$pos - Start)) %>%
      filter(dis < dis_cis/2) %>%
      arrange(dis)
    
    c(
      "nearest_dis_module" = tmp_coexp_module$dis[1],
      "nearest_gene_module" = tmp_coexp_module$gene[1],
      "cis_dis_module" = paste(tmp_coexp_module$dis, collapse = ";"),
      "cis_gene_module" = paste(tmp_coexp_module$gene, collapse = ";")
    )
  }
)

cis_gene_module <- as.data.table(t(cis_gene_module)) %>%
  mutate(nearest_dis_module = as.numeric(nearest_dis_module))


signal_cis_genes <- bind_cols(signal_cis_genes, cis_gene_module) %>%
  relocate(nearest_dis_module, nearest_gene_module, .after = p) %>%
  relocate(cis_dis_module, cis_gene_module, .after = nearest_gene_module) %>%
  arrange(nearest_dis_module)

fwrite(signal_cis_genes, 'dis_to_genes_in_module.txt', quote = FALSE, sep = "\t")
