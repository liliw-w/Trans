###############################################################
########### Prepare gene set across module ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
file_gene_set <- '/project2/xuanyao/llw/ldsc/test_ldsc/test.GeneSet'

# output -----
file_gene_ensg <- '/project2/xuanyao/llw/ldsc/test_ldsc/test_ensg.GeneSet'


# read data -----
gene_meta <- fread(file_gene_meta, header = TRUE)
gene_set <- fread(file_gene_set, header = FALSE, col.names = "gene")


# 1. extract genes for each module -----
gene_set <- left_join(gene_set, gene_meta, by = c("gene" = "GeneSymbol")) %>%
  separate(Geneid, c("ensg", NA), "[.]")


# 2. keep only trans genes in the module for a chr
#gene_trans <- gene_meta %>%
#  filter(!is.na(module) & module == !!module & chr != paste0("chr", !!chr)) %>%
#  select(ensg)


# output: geneset for each module -----
fwrite(
  gene_set$ensg %>% enframe(name = NULL),
  file = file_gene_ensg,
  quote = FALSE, sep = "\t", col.names = FALSE
)

