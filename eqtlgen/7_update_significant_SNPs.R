###############################################################
########### Update eQTLGen signals with the new BF correction ###########
########### which uses only the modules after ratio thresholding ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


### paras and I/O
ratio <- 50
thre_p_z <- 1e-4
n_SNP <- 9918
fdr_level <- 0.1


file_p_all <- 'p/p.module_all.Sigma_nullz.rds'
file_module_use <- paste0('postanalysis/module_use_ratio_', ratio, '.txt')

file_signal <- paste0('postanalysis/signal_rm_infl_ratio_', ratio, '.txt')
file_signal_uniq <- paste0('postanalysis/signal_uniq_rm_infl_ratio_', ratio, '.txt')
file_signal_LD <- 'postanalysis/LD.prun.in.chr.module.txt'


### read files
p_all <- readRDS(file_p_all)
module_use <- fread(file_module_use, header = TRUE)


### filter out module with #NULL SNPs/Module size ratio less than the given ratio
module_use <- pull(module_use, module)


### use Bonforroni correction based on the remained modules
n_module <- length(module_use)
p_thre <- fdr_level/n_SNP/n_module


### extract signals
signal <- p_all %>%
  filter(p <= p_thre & module %in% module_use) %>%
  arrange(SNPChr, SNPPos, p)
signal_uniq <- signal %>% distinct(SNP, .keep_all = TRUE) %>%
  select(SNP, SNPChr, SNPPos, meta)


### write out signals
fwrite(signal,
       file_signal,
       sep = "\t", quote = FALSE)

fwrite(signal_uniq,
       file_signal_uniq,
       sep = "\t", quote = FALSE)

fwrite(signal_uniq %>% distinct(meta),
       file_signal_LD,
       sep = "\t", quote = FALSE, col.names = FALSE)


