##############################################
########### Make g:profiler enrichment file for multiple inqueries ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_coexp_net <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds'

## output -----
dir_module_enrich_query <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/'

# read files -----
coexp_net <- readRDS(file_coexp_net)


# write out in format of g:profile multiple query -----
enframe(coexp_net$moduleLabels, name = "gene", value = "module") %>%
  group_by(module) %>%
  summarise(
    gene = paste(gene, collapse = " ")
  ) %>%
  ungroup() %>%
  filter(module != 0) %>%
  mutate(
    g = cut(module, c(seq(1, max(module), by = 20), Inf), right = FALSE),
    query = paste0("> M", module),
    module = NULL
  ) %>%
  pivot_longer(c(query, gene), names_to = NULL, values_to = "enrich") %>%
  group_by(g) %>%
  group_walk(~ fwrite(
    .x, paste0(dir_module_enrich_query, "query_group_", .y$g, ".txt"), 
    sep = '\t', quote = FALSE, col.names = FALSE
  ) ) %>%
  ungroup()

