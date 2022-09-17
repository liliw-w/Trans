###############################################################
########### claim significant SNPs ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_p_all <- 'FDR/p_all_1e5.rds'
p_thre <- 0.1/9918/50

## output -----
file_signal <- 'postanalysis/signal.txt'
file_signal_uniq <- 'postanalysis/signal_uniq.txt'


# read files -----
p_all <- readRDS(file_p_all)


# extract signals -----
signal <- filter(p_all, p <= p_thre) %>%
  separate(module, c(NA, "module"), sep = "module", convert = TRUE) %>%
  select(module, meta, p, SNP, SNPChr, SNPPos) %>%
  arrange(module, SNPChr, SNPPos, p)
signal_uniq <- distinct(signal, SNP, .keep_all = TRUE) %>%
  select(SNP, SNPChr, SNPPos, meta)


# print out key message or write out -----
fwrite(signal,
       file_signal,
       sep = "\t", quote = FALSE)

fwrite(signal_uniq,
       file_signal_uniq,
       sep = "\t", quote = FALSE)
