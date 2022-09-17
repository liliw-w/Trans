###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Remove modules with ratio under a specified value ###########
########### Write out modules used for follow up analysis without signal inflation ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


### paras and I/O
thre_p_z <- 1e-4
ratio <- 100

file_null_SNP <- 'null_SNP/num_nullSNP.rds'

file_module_use <- paste0('postanalysis/module_use_ratio_', ratio, '.txt')


### read files
res_nullSNP <- readRDS(file_null_SNP)


### select modules with ratio higher than the given ratio
module_use <- res_nullSNP %>%
  filter(thre_z == thre_p_z) %>%
  mutate(prop_dim = num_nullSNP_indep/module_size) %>%
  filter(prop_dim > ratio)


fwrite(module_use,
       file_module_use,
       sep = "\t", quote = FALSE)
