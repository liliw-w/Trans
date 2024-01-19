rm(list = ls())
library(tidyverse)


# 0. prep data for plotting -----
## read data -----
file_target_enrich <- list.files(".", "*target_enrich_hallmark.txt", full.names = TRUE)

target_enrich <- lapply(
  file_target_enrich,
  FUN = data.table::fread
) %>%
  bind_rows()



## use genes for plotting -----
plt_dat_gene <- select(
  target_enrich,
  trait, nearest_gene, if_trans_sig, if_target_all, if_target_trait, if_immune_blood_anno
) %>%
  group_by(trait, nearest_gene) %>%
  summarise(
    if_trans_sig = any(if_trans_sig),
    if_target_all = any(if_target_all),
    if_target_trait = any(if_target_trait),
    if_immune_blood_anno = any(if_immune_blood_anno)
  ) %>%
  ungroup() %>%
  count(trait, if_trans_sig, if_target_all, if_target_trait, if_immune_blood_anno)

plt_dat_gene$if_trans_sig_immune_blood_anno <- plt_dat_gene$if_trans_sig & as.logical(plt_dat_gene$if_immune_blood_anno)



dat_enrich <- group_by(plt_dat_gene, trait, if_target_all, if_trans_sig_immune_blood_anno) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  arrange(
    trait, desc(if_target_all), desc(if_trans_sig_immune_blood_anno)
  )

res_fisher <- list()
for(k in unique(dat_enrich$trait)){
  tmp_dat_enrich = filter(dat_enrich, trait == !!k)
  
  if(nrow(tmp_dat_enrich) == 4) {
    tmp_fisher = fisher.test(matrix(tmp_dat_enrich$n, ncol = 2), alternative = "greater")
    
    res_fisher = c(
      res_fisher,
      list(tibble(
        'trait' = k,
        'p_enrich' = tmp_fisher$p.value, 
        'or_enrich' = tmp_fisher$estimate
      ))
      )
  }
}
res_fisher <- bind_rows(res_fisher)

res_fisher

