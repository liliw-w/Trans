##############################################
########### look into what are those pc1 signals missed by pco ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
file_pco_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_pc1_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability_PC1/FDR/signals.chr.module.perm10.txt'
file_pco_p_all <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/q.chr.module.perm1.rds'

file_mis_sig <- "pc1/pc1_sig_mis_pco.txt"


# read files -----
pco_sig <- fread(file_pco_sig, header = FALSE, col.names = c('signal', 'p', 'q'))
pc1_sig <- fread(file_pc1_sig, header = FALSE, col.names = c('signal', 'p', 'q'))
pco_p_all <- readRDS(file_pco_p_all)


# which are pc1 signals missed by pco & add p and q for pc1 and pco -----
pc1_sig <- pc1_sig %>%
  mutate(
    "if_pc1" = TRUE,
    "if_pco" = signal %in% !!pco_sig$signal
  ) %>%
  left_join(pco_p_all, by = c("signal" = "snp"), , suffix = c("_pc1", "_pco")) %>%
  relocate(starts_with("if_"), .after = signal) %>%
  arrange(p_pc1)


fwrite(
  pc1_sig %>% filter(if_pc1 & !if_pco),
  file_mis_sig,
  quote = FALSE, sep = "\t"
)
