###############################################################
########### claim significant SNPs ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)

file_p_all <- 'p/p.module_all.Sigma_nullz.rds'
file_signal <- 'postanalysis/signal.txt'
file_signal_uniq <- 'postanalysis/signal_uniq.txt'

### use Bonforroni correction
p_thre <- 0.1/9918/166


### read data
p_all <- readRDS(file_p_all)


### extract signals
signal <- p_all %>% filter(p <= p_thre)
signal_uniq <- signal %>% distinct(SNP, .keep_all = TRUE) %>%
  select(SNP, SNPChr, SNPPos, meta)


### write out signals
fwrite(signal,
       file_signal,
       sep = "\t", quote = FALSE)

fwrite(signal_uniq,
       file_signal_uniq,
       sep = "\t", quote = FALSE)
