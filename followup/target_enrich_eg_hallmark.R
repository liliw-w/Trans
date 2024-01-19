##############################################
########### Plot to check if drug targets are more likely to be trans signals of gene modules ###########
########### or if trans signals of gene modules are more likely to be drug targets ###########
########### use hallmark gene sets as modules ###########
##############################################
rm(list = ls())
library(tidyverse)



# read data -----
file_target_enrich <- list.files(".", "*target_enrich_hallmark.txt", full.names = TRUE)
file_module_meta_enrich <- '/project2/xuanyao/llw/MODULES/MSigDB/msig_50_table_immune_related.txt'


target_enrich <- lapply(
  file_target_enrich,
  FUN = data.table::fread
) %>%
  bind_rows()
module_meta_enrich <- data.table::fread(file_module_meta_enrich)



# 1. targets genes that are trans -----
filter(target_enrich, if_target_all & if_trans_sig) %>%
  distinct(nearest_gene) %>%
  arrange(nearest_gene)

filter(target_enrich, if_target_all & if_trans_sig) %>%
  select(
    nearest_gene, if_immune_blood_anno, trait, gene_module_immune
  )


# 2. targets genes that are trans corresponding to immune-related modules -----
filter(target_enrich, if_target_all & if_trans_sig & if_immune_blood_anno) %>%
  distinct(nearest_gene) %>%
  arrange(nearest_gene)



# 3. trans module annotation -----
module_of_interest <- filter(target_enrich, if_target_all & if_trans_sig) %>%
  pull(gene_module_immune) %>%
  str_split(";") %>%
  unlist() %>%
  unique()
module_of_interest


filter(module_meta_enrich, TERMS %in% module_of_interest) %>%
  select(if_immune_blood_anno, TERMS, description_info)

