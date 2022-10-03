###############################################################
########### Update eQTLGen signals ###########
########### by removing modules and chr's with ratio lower than cutoff ###########
###############################################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
ratio_thre <- 50
thre_p_z <- 1e-4
fdr_level <- 0.1
n_snp <- 9918


file_p_all <- 'FDR/p_all_1e5.rds'
file_ratio_module_chr <- 'null_SNP/ratio_module_chr.txt'

## output -----
file_signal_rm_infl_ratio <- paste0('postanalysis/signal_rm_infl_ratio_', ratio_thre, '.txt')
file_signal_uniq_rm_infl_ratio <- paste0('postanalysis/signal_uniq_rm_infl_ratio_', ratio_thre, '.txt')


# read files -----
p_all <- readRDS(file_p_all)
ratio_module_chr <- fread(file_ratio_module_chr, header = TRUE)


# filter out modules and chr's with #NULL SNPs/Module size ratio less than the given ratio -----
module_chr_use <- filter(ratio_module_chr, prop_dim_module > ratio_thre) %>%
  distinct(module, SNPChr, prop_dim_module)


# determin p cutoff based on chosen modules -----
p_thre <- fdr_level/n_snp/n_distinct(module_chr_use$module)


# extract signals based on ratio cutoff -----
# extract signals -----
signal <- filter(p_all, p <= p_thre) %>%
  separate(module, c(NA, "module"), sep = "module", convert = TRUE) %>%
  select(module, meta, p, SNP, SNPChr, SNPPos) %>%
  arrange(module, SNPChr, SNPPos, p) %>%
  inner_join(
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

