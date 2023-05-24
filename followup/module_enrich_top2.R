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
  pattern = '^enrich_query_group_.*[.]csv$'
)

## output -----
file_enrich_top2 <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich_top2.txt'



# extract the top2 enrichment of term categories in modules -----
enrich_top2 <- NULL
for(file_enrich in file_enrich_list){
  dat <- fread(file_enrich, sep = ",") %>% select(!p_adj_color)
  nmod <- (ncol(dat) - 5) / 3
  
  for(i in 1:nmod){
    tmp_dat <- select(dat, c(1, 2, i*3 + 3))
    M <- str_extract(colnames(tmp_dat)[ncol(tmp_dat)], '\\d+$')
    
    tmp_enrich_top2 <- rename(tmp_dat, adjusted_p_value = last_col()) %>%
      filter(adjusted_p_value < 0.05 & adjusted_p_value >=0) %>%
      group_by(source) %>%
      arrange(adjusted_p_value, .by_group = TRUE) %>%
      slice(1:2) %>%
      ungroup() %>%
      pivot_wider(
        names_from = source, 
        values_from = !source, 
        names_glue = "{source}_{.value}", 
        values_fn = ~paste(.x, collapse = ";")
      ) %>%
      mutate(gene_module = M) %>%
      relocate(gene_module, .before = everything())
    enrich_top2 <- bind_rows(enrich_top2, tmp_enrich_top2)
  }
}


# print out key message or write out -----
fwrite(
  enrich_top2,
  file_enrich_top2,
  sep = '\t', quote = FALSE
)

