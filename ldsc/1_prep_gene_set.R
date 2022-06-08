###############################################################
########### Prepare gene set across module ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
file_gene_meta <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/gene.meta.txt'
file_coexp_module <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds'

dir_out <- 'geneset'


# read data -----
gene_meta <- fread(file_gene_meta, header = TRUE)
coexp_module <- readRDS(file_coexp_module)$moduleLabels


# 1. extract genes for each module -----
gene_meta <- enframe(coexp_module, "gene", "module") %>%
  right_join(gene_meta, by = c("gene")) %>%
  separate(GeneNameConv, c("ensg", NA), "[.]")

# 2. keep only trans genes in the module for a chr
#gene_trans <- gene_meta %>%
#  filter(!is.na(module) & module == !!module & chr != paste0("chr", !!chr)) %>%
#  select(ensg)


# output: geneset for each module -----
gene_meta %>%
  filter(!is.na(module) & module != 0) %>%
  group_by(module) %>%
  group_walk(
    ~fwrite(.x$ensg %>% enframe(name = NULL),
            file = file.path(dir_out, paste0("M", .y$module, ".GeneSet")),
            quote = FALSE, sep = "\t", col.names = FALSE)
  )
