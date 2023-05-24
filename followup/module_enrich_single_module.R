##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_enrich_list <- list.files(
  '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich',
  pattern = '^enrich_query_group_.*[.]csv$',
  full.names = TRUE
)

## output -----
dir_enrich_single_M <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/enrich_single_module/'



# extract the single_M enrichment of term categories in modules -----
enrich_single_M <- NULL
for(file_enrich in file_enrich_list){
  dat <- fread(file_enrich, sep = ",") %>% select(!p_adj_color)
  nmod <- (ncol(dat) - 5) / 3
  
  for(i in 1:nmod){
    tmp_dat <- select(dat, c(1, 2, i*3 + 3))
    M <- str_extract(colnames(tmp_dat)[ncol(tmp_dat)], '\\d+$')
    
    tmp_enrich_single_M <- rename(tmp_dat, adjusted_p_value = last_col()) %>%
      filter(adjusted_p_value < 0.05 & adjusted_p_value >=0) %>%
      group_by(source) %>%
      arrange(adjusted_p_value, .by_group = TRUE) %>%
      ungroup() %>%
      mutate(gene_module = M) %>%
      relocate(gene_module, .before = everything())
    enrich_single_M <- bind_rows(enrich_single_M, tmp_enrich_single_M)
  }
}


# save module enrichment for every single module -----
enrich_single_M %>%
  group_by(gene_module) %>%
  group_walk(~ fwrite(
    .x, paste0(dir_enrich_single_M, "module_enrich_M", .y$gene_module, ".txt"), 
    sep = '\t', quote = FALSE
  ), .keep = TRUE ) %>%
  ungroup()

