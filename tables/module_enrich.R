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
file_enrich_all <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich.txt'


# extract all enrichment of term categories in modules -----
enrich_all <- NULL
for(file_enrich in file_enrich_list){
  dat <- fread(file_enrich, sep = ",") %>% select(!p_adj_color)
  nmod <- (ncol(dat) - 5) / 3
  
  for(i in 1:nmod){
    tmp_dat <- select(dat, c(1, 2, 3, i*3 + 3))
    M <- str_extract(colnames(tmp_dat)[ncol(tmp_dat)], '\\d+$')
    
    tmp_enrich_all <- rename(tmp_dat, adjusted_p_value = last_col()) %>%
      filter(adjusted_p_value < 0.05 & adjusted_p_value >=0) %>%
      group_by(source) %>%
      arrange(adjusted_p_value, .by_group = TRUE) %>%
      ungroup() %>%
      pivot_wider(
        names_from = source, 
        values_from = !source, 
        names_glue = "{source}_{.value}", 
        values_fn = ~paste(.x, collapse = ";")
      ) %>%
      mutate(gene_module = M) %>%
      relocate(gene_module, .before = everything())
    
    enrich_all <- bind_rows(enrich_all, tmp_enrich_all)
  }
}


# adjust row & col order -----
enrich_all <- select(
  enrich_all,
  gene_module, 
  `GO:BP_term_name`, `GO:BP_term_id`, `GO:BP_adjusted_p_value`, 
  `GO:MF_term_name`, `GO:MF_term_id`, `GO:MF_adjusted_p_value`, 
  `KEGG_term_name`, `KEGG_term_id`, `KEGG_adjusted_p_value`, 
  `REAC_term_name`, `REAC_term_id`, `REAC_adjusted_p_value`
) %>%
  mutate(gene_module = as.numeric(gene_module)) %>%
  arrange(gene_module)


# print out key message or write out -----
fwrite(
  enrich_all,
  file_enrich_all,
  sep = '\t', quote = FALSE
)

