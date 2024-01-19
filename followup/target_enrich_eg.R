##############################################
########### Plot to check if drug targets are more likely to be trans signals of gene modules ###########
########### or if trans signals of gene modules are more likely to be drug targets ###########
##############################################
rm(list = ls())
library(tidyverse)

trait_of_interest <- "allergy"


# read data -----
file_target_enrich <- list.files(".", "*target_enrich.txt", full.names = TRUE)
file_module_meta_enrich <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich_top2_immune_related.txt'


target_enrich <- lapply(
  file_target_enrich,
  FUN = data.table::fread
) %>%
  bind_rows()
module_meta_enrich <- data.table::fread(file_module_meta_enrich)



# 1. targets genes that are trans -----
target_enrich <- filter(
  target_enrich,
  trait == !!trait_of_interest
)


filter(target_enrich, if_target_all & if_trans_sig) %>%
  distinct(nearest_gene) %>%
  arrange(nearest_gene)



# 2. targets genes that are trans corresponding to immune-related modules -----
filter(target_enrich, if_target_all & if_trans_sig & if_immune_blood_anno) %>%
  distinct(nearest_gene) %>%
  arrange(nearest_gene)



# 3. trans module annotation -----
module_of_interest <- filter(target_enrich, if_target_all & if_trans_sig) %>%
  pull(gene_module_immune) %>%
  str_split(";") %>%
  unlist() %>%
  unique() %>%
  as.numeric()

module_of_interest


filter(module_meta_enrich, gene_module %in% module_of_interest) %>%
  select(gene_module, ends_with("_term_name"))

