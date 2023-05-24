##############################################
########### numeric results of coloc prop for traits ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_res_coloc_reg_prop1 <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/coloc_region_prop_merged.txt'
file_res_coloc_reg_prop2 <- '/project2/xuanyao/llw/coloc/immune_traits/pmid_all/coloc_region_prop_merged.txt'
file_res_coloc_reg_prop3 <- '/project2/xuanyao/llw/coloc/ukbb_coloc_more_traits/all_trait/data/coloc_region_prop_merged.txt'

file_blood_meta <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv'


## output -----
file_num_prop <- "/project2/xuanyao/llw/coloc/coloc_prop_all.txt"


# read files -----
res_coloc_reg_prop1 <- fread(file_res_coloc_reg_prop1, header = TRUE)
res_coloc_reg_prop2 <- fread(file_res_coloc_reg_prop2, header = TRUE)
res_coloc_reg_prop3 <- fread(file_res_coloc_reg_prop3, header = TRUE)
blood_meta <- fread(file_blood_meta) %>% mutate(`GWAS ID` = as.character(`GWAS ID`))


# change order and rename cols -----
res_coloc_reg_prop1 <- mutate(res_coloc_reg_prop1, trait_type = "Blood") %>% arrange(desc(propPvalColocMerg))
res_coloc_reg_prop2 <- mutate(res_coloc_reg_prop2, trait_type = "Autoimmune") %>% arrange(desc(propPvalColocMerg))
res_coloc_reg_prop3 <- mutate(res_coloc_reg_prop3, trait_type = "Other") %>% arrange(desc(propPvalColocMerg))

res_coloc_reg_prop <- rbind(res_coloc_reg_prop1, res_coloc_reg_prop2, res_coloc_reg_prop3) %>%
  left_join(select(blood_meta, `GWAS ID`, `GWAS Trait`), by = c('Phenocode' = 'GWAS ID')) %>%
  select(
    trait_type, Phenocode, trait, `GWAS Trait`, 
    nRegionPval, nRegionPvalColoc, propPvalColoc, 
    nRegionPvalMerg, nRegionPvalColocMerg, propPvalColocMerg
  ) %>%
  rename(
    phenocode = Phenocode, trait_info = `GWAS Trait`,
    n_region = nRegionPval, n_region_coloc = nRegionPvalColoc, prop_coloc = propPvalColoc, 
    n_region_merg = nRegionPvalMerg, n_region_merg_coloc = nRegionPvalColocMerg, prop_merg_coloc = propPvalColocMerg
  )

# print out key message or write out -----
fwrite(
  res_coloc_reg_prop,
  file_num_prop,
  sep = '\t', quote = FALSE
)

