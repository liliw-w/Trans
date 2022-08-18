rm(list = ls())
library(data.table)
library(tidyverse)

pp4Thre = 0.75
pvalThre = 'module_QTL_sig'
nsnpsThre = 5

dis_cis_diff <- 2e+5

dir_coloc = file.path( "~/xuanyao_llw/coloc/ukbb_coloc_blood_traits")


file_resColoc = paste0("/scratch/midway2/liliw1/coloc/data/pheno", gwasPhenocode, ".resColoc.txt.gz")


file_resColoc = list.files(dir_coloc, ".*resColoc.txt.gz", full.names = TRUE, recursive = TRUE)
file_module_QTL_signals = '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'

#file_resColoc <- file_resColoc_list[1]

module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))

resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)
resColoc_all = resColoc_all %>% mutate('if_module_QTL' = Region %in% module_QTL_signals$signal )


res_coloc_reg_prop = resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegionPval = sum(if_module_QTL),
            nRegionPvalColoc = sum(if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)) %>%
  ungroup()


res_coloc_reg_merged <- resColoc_all %>%
  filter(if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) %>%
  separate(Region, into = c("M", "chr", "pos"), sep = "[:]", remove = FALSE, convert = TRUE) %>%
  group_by(Phenocode, trait, M, chr) %>%
  summarise(
    "n_reg" = n(),
    "dis" = paste(diff(sort(pos)), collapse = ";"),
    "if_cis" = paste(diff(sort(pos)) < dis_cis_diff, collapse = ";"),
    "num_merged_reg" = sum(diff(sort(pos)) > dis_cis_diff) + 1,
    "comb_reg" = paste(Region, collapse = ";")
  ) %>%
  ungroup() %>%
  group_by(Phenocode, trait) %>%
  summarise(
    nRegionPvalColoc_merged = sum(num_merged_reg)
  ) %>%
  ungroup() %>%
  right_join(
    res_coloc_reg_prop,
    by = c("Phenocode", "trait")
  ) %>%
  relocate(
    nRegionPvalColoc_merged,
    .after = everything()
  ) %>%
  arrange(
    desc(nRegionPvalColoc_merged)
  )

fwrite(res_coloc_reg_merged, "res_coloc_reg_merged.txt", sep = "\t")


res_coloc_reg_prop = res_coloc_reg_prop %>%
  mutate(propColoc = nRegionColoc/nRegion,
         propPvalColoc = nRegionPvalColoc/nRegionPval) %>%
  arrange(desc(propPvalColoc), desc(propColoc)) 


# file 1: coloc & if_module_QTL across each trait



# file 2: merged coloc regions across each trait

# file 3: coloc & if_module_QTL for all traits

# file 4: merged coloc regions for all traits

