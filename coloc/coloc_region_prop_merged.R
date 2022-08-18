############################################################
########## Merged coloc regions & proportions for traits ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_list_resColoc <- list.files(
  "/scratch/midway2/liliw1/tmp",
  "^pheno\\d+[.]coloc_reg_w_merged.txt$",
  full.names = TRUE
)
file_reg_merged <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/qtl_region_merged.txt"

file_res_coloc_reg_prop <- 'coloc_region_prop_merged.txt'


# read files -----
resColoc_all <- rbindlist(lapply(file_list_resColoc, fread, header = TRUE), use.names = TRUE)
reg_merged <- fread(file_reg_merged, header = TRUE)


# add coloc regions, unmerged and merged. also add coloc proportion -----
res_coloc_reg_prop <- resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegionPval = n_distinct(!!reg_merged$Region),
            nRegionPvalColoc = n_distinct(Region),
            nRegionPvalMerg = n_distinct(!!reg_merged$merged_region),
            nRegionPvalColocMerg = n_distinct(merged_region)) %>%
  ungroup() %>%
  mutate(propPvalColoc = nRegionPvalColoc/nRegionPval,
         propPvalColocMerg = nRegionPvalColocMerg/nRegionPvalMerg) %>%
  arrange(desc(propPvalColocMerg))


fwrite(
  res_coloc_reg_prop,
  file = file_res_coloc_reg_prop,
  quote = FALSE,
  sep = '\t'
)
