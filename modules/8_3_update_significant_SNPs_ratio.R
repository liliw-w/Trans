###############################################################
########### Update eQTLGen signals ###########
########### by removing modules and chr's with ratio lower than cutoff ###########
###############################################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
thre_p_z <- 1e-4
ratio_thre <- 50

file_signal <- 'postanalysis/signal.txt'
file_ratio_module_chr <- 'null_SNP/ratio_module_chr.txt'

## output -----
file_signal_rm_infl_ratio <- paste0('postanalysis/signal_rm_infl_ratio_', ratio_thre, '.txt')
file_signal_uniq_rm_infl_ratio <- paste0('postanalysis/signal_uniq_rm_infl_ratio_', ratio_thre, '.txt')


# read files -----
signal <- fread(file_signal, header = TRUE)
ratio_module_chr <- fread(file_ratio_module_chr, header = TRUE)


# filter out modules and chr's with #NULL SNPs/Module size ratio less than the given ratio -----
module_chr_use <- filter(ratio_module_chr, prop_dim_module > ratio_thre) %>%
  distinct(module, SNPChr, prop_dim_module)


# extract signals based on ratio cutoff -----
signal <- inner_join(
  signal,
  module_chr_use,
  by = c("module", "SNPChr")
) %>%
  select(module, meta, p, SNP, SNPChr, SNPPos) %>%
  arrange(module, SNPChr, SNPPos, p)

signal_uniq <- distinct(signal, SNP, .keep_all = TRUE) %>%
  select(SNP, SNPChr, SNPPos, meta)


# print out key message or write out -----
fwrite(signal,
       file_signal_rm_infl_ratio,
       sep = "\t", quote = FALSE)

fwrite(signal_uniq,
       file_signal_uniq_rm_infl_ratio,
       sep = "\t", quote = FALSE)

