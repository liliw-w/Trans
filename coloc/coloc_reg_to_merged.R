############################################################
########## Assign coloc regions of a trait to merged reion ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)

# I/O & paras -----
pp4Thre <- 0.75
pvalThre <- 'module_QTL_sig'
nsnpsThre <- 5

file_list_resColoc <- list.files(
  "~/xuanyao_llw/coloc/ukbb_coloc_blood_traits/data",
  "^pheno\\d+[.]resColoc.txt.gz",
  full.names = TRUE
)
file_module_QTL_signals <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_reg_merged <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/qtl_region_merged.txt"


# read files -----
module_QTL_signals <- fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
reg_merged <- fread(file_reg_merged, header = TRUE)


# select coloc region and assign it to a merged region -----
for(file_resColoc in file_list_resColoc){
  cat("Running", file_resColoc, '\n\n')
  
  resColoc <- fread(file_resColoc, header = TRUE)
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
    file = str_replace(basename(file_resColoc), "resColoc.txt.gz$", "coloc_reg_w_merged.txt"),
    quote = FALSE,
    sep = '\t'
  )
}
