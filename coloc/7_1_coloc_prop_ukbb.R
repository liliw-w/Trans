##############################################
########### calcualte coloc proportion ##########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
pp4Thre <- 0.75
pvalThre <- 'module_QTL_sig'
nsnpsThre <- 5

file_module_QTL_signals <- '/project2/xuanyao/llw/MODULES/MSigDB/FDR/signals.chr.module.perm10.fdr10.txt'
dir_coloc_gwas <- list.files(
  "/scratch/midway2/liliw1/coloc_MSigDB",
  pattern = paste0("^ukbb_continuous_\\d+"),
  full.names = TRUE
)
file_resColoc <- sapply(
  dir_coloc_gwas,
  list.files,
  pattern = "^resColoc.txt.gz$", full.names = TRUE, recursive = TRUE
)


## output -----
file_res_coloc_reg_prop <- paste0("ukbb_all/coloc_region_prop_pvalThre-", pvalThre, ".txt")
file_coloc_all <- paste0("ukbb_all/coloc_region_summary_all_pvalThre-", pvalThre, ".txt")


# read files -----
module_QTL_signals <- fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
range(module_QTL_signals$p)
range(module_QTL_signals$q)

resColoc_all <- rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)


# number of regions & coloc regions, proportion of coloc regions -----
resColoc_all <- mutate(
  resColoc_all,
  'if_module_QTL' = Region %in% module_QTL_signals$signal
)

res_coloc_reg_prop <- resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegion = n(),
            nRegionColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre),
            nRegionPval = sum(if_module_QTL),
            nRegionPvalColoc = sum(if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)) %>%
  ungroup()
res_coloc_reg_prop <- res_coloc_reg_prop %>%
  mutate(propColoc = nRegionColoc/nRegion,
         propPvalColoc = nRegionPvalColoc/nRegionPval) %>%
  arrange(desc(propPvalColoc), desc(propColoc))


# print out key message or write out -----
fwrite(
  filter(resColoc_all, if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre),
  file_coloc_all,
  quote = FALSE, sep = "\t"
)
fwrite(res_coloc_reg_prop, file_res_coloc_reg_prop, quote = FALSE, sep = "\t")

