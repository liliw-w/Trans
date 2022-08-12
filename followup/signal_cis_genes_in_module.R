rm(list = ls())
library(data.table)
library(tidyverse)

file_signal_cis_genes = '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/signal_cis_genes.txt'
file_coexp_module <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds'

file_cis_genes_in_module <- 'cis_genes_in_module.txt'


signal_cis_genes <- fread(file_signal_cis_genes, header = TRUE)
coexp_module <- readRDS(file_coexp_module)$moduleLabels

coexp_module <- enframe(coexp_module, "gene", "m")


cis_gene_module <- sapply(
  1:nrow(signal_cis_genes),
  function(x){
    tmp_cis = signal_cis_genes[x, ]
    tmp_coexp_module = coexp_module %>% filter(paste0("module", m) == tmp_cis$module)
    
    nearest_gene = tmp_cis$nearest_gene
    near_genes = tmp_cis$near_genes %>% str_split(pattern = ";") %>% unlist()
    
    c(
      sum(nearest_gene %in% tmp_coexp_module$gene),
      paste0("M", coexp_module %>% filter(gene == nearest_gene) %>% pull(m)),
      
      sum(near_genes %in% tmp_coexp_module$gene),
      enframe(near_genes) %>%
        left_join(coexp_module, by = c("value" = "gene")) %>%
        mutate(m = as.character(m)) %>%
        replace_na(list(m = "")) %>%
        pull(m) %>%
        paste0("M", ., collapse = ";")
    )
  }
)
cis_gene_module <- as.data.table(t(cis_gene_module))
colnames(cis_gene_module) = c("if_in_module", "module_gene", "if_in_module_near", "module_genee_near")


signal_cis_genes <- bind_cols(signal_cis_genes, cis_gene_module) %>%
  relocate(if_in_module, module_gene, .after = nearest_gene) %>%
  relocate(if_in_module_near, module_genee_near, .after = near_genes)

fwrite(signal_cis_genes, file_cis_genes_in_module, quote = FALSE, sep = "\t")
