############################################################
########## Assign coloc regions of a trait to merged region ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)

# I/O & paras -----
pp4Thre <- 0.75
pvalThre <- 'module_QTL_sig'
nsnpsThre <- 5

file_resColoc <- "/project2/xuanyao/llw/coloc/cis/cis_e/data/resColoc.txt.gz"
file_module_QTL_signals <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_reg_merged <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/qtl_region_merged.txt"

## output -----
file_resColoc_merged <- "/project2/xuanyao/llw/coloc/cis/cis_e/data/coloc_reg_w_merged.txt"


# read files -----
resColoc <- fread(file_resColoc, header = TRUE)
module_QTL_signals <- fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
reg_merged <- fread(file_reg_merged, header = TRUE)


# select coloc region and assign it to a merged region -----
resColoc <- resColoc %>%
  filter(Region %in% module_QTL_signals$signal &
           PP.H4.abf > pp4Thre &
           nsnps >= nsnpsThre) %>%
  left_join(reg_merged, by = "Region") %>%
  arrange(merged_region) %>%
  relocate(
    merged_region, merged_region_n, PP.H4.abf,
    .after= Region
  )

fwrite(
  resColoc,
  file = file_resColoc_merged,
  quote = FALSE,
  sep = '\t'
)
