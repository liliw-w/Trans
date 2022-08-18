############################################################
########## Regions for coloc defined by trans-eQTLs ##########
########## and merged regions based on that ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)

# I/O & paras -----
file_qtl_reg <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/qtlColocReg.txt.gz'
file_module_QTL_signals <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'

dis_cis_diff <- 2e+5

file_reg_merged <- "qtl_region_merged.txt"


# read files -----
qtl_reg <- fread(file_qtl_reg)
module_QTL_signals <- fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))


# merge regions with close lead SNPs for coloc -----
reg_merged <- qtl_reg %>%
  filter(
    Region %in% !!module_QTL_signals$signal,
    Signal == Region
  ) %>%
  group_by(Module, Chr) %>%
  arrange(Pos, .by_group = TRUE) %>%
  group_modify(~ {
    cut(seq(1, length(.x$Pos)),
        breaks = c(which(c(TRUE, diff(.x$Pos) > dis_cis_diff)), Inf),
        right = FALSE) %>%
      enframe(name = NULL, value = "l") %>%
      cbind('Region' = .x$Region, "Pos" = .x$Pos) %>%
      group_by(l) %>%
      mutate("merged_region" = paste(Region, collapse = ";"),
             "merged_region_n" = n())
  }) %>%
  ungroup() %>%
  select(!l)


# write out -----
fwrite(
  reg_merged,
  file_reg_merged,
  quote = FALSE,
  sep = "\t"
)

cat(
  "There were", nrow(reg_merged), "regions. \n\n",
  "After merging close lead SNPs of the regions, there are",
  nrow(distinct(reg_merged, merged_region)), "merged regions. \n\n"
)
