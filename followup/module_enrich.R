###############################################################
########### bar plot coloc with cis-e/s numbers ###########
########### use only numbers appear in manuscript ###########
###############################################################
rm(list = ls())
library(tidyverse)


# paras and I/O -----
file_coexp_net <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/coexp.module.rds'


# read coexpression network -----
coexp_net <- readRDS(file_coexp_net)


# write module genes into separate files -----
M_df <- enframe(coexp_net$moduleLabels, name = "gene", value = "M")

count(M_df, M) %>% arrange(M)

M_df %>% group_by(M) %>% group_walk(
  ~write_tsv(paste(.x$gene, collapse = " ") %>% tibble(), paste0("module_enrich/module", .y, ".txt"), col_names = FALSE)
  )

