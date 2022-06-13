###############################################################
###########  Look at the replication between DGN signals and eQTLGen trans signals ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


### paras and I/O
fdr_level <- 0.1
ratio <- 50


file_dgn_all_snp <- "/project2/xuanyao/llw/eQTLGen_DGN/DGN.all_snp.txt"
file_eqtlgen_all_snp <- "/project2/xuanyao/llw/eQTLGen/meta.snp.txt.gz"

file_dgn_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_eqtlgen_old_signal <- paste0('postanalysis/signal_rm_infl_ratio_', ratio, '.txt')

file_p_all <- 'p/p.module_all.Sigma_nullz.rds'
file_module_use <- paste0('postanalysis/module_use_ratio_', ratio, '.txt')


### read files
dgn_all_snp <- fread(file_dgn_all_snp, header = FALSE)
dgn_sig <- fread(file_dgn_sig, header = FALSE, col.names = c("signal", "p", "q") )

eqtlgen_all_snp <- fread(file_eqtlgen_all_snp, header = TRUE) # 10317 all eQTLGen SNPs
eqtlgen_old_signal <- fread(file_eqtlgen_old_signal, header = TRUE)

p_all <- readRDS(file_p_all)
module_use <- fread(file_module_use, header = TRUE)


### extract SNPs
dgn_all_snp_uniq <- dgn_all_snp$V1

eqtlgen_all_snp <- unite(eqtlgen_all_snp, "SNP_id", c(SNPChr, SNPPos), sep = ":", remove = FALSE)
eqtlgen_all_snp_uniq <- eqtlgen_all_snp$SNP_id


dgn_sig <- dgn_sig %>%
  separate(signal,
           into = c('module', 'chr', 'pos'),
           sep = ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_id", c(chr, pos), sep = ":", remove = FALSE)
dgn_sig_uniq <- dgn_sig %>% distinct(SNP_id) %>% pull()


n_dgn_sig_rm <- dgn_sig %>%
  filter(module %in% paste0("module", module_use$module) &
           SNP_id %in% eqtlgen_all_snp$SNP_id ) %>%
  nrow()

p_thre <- fdr_level/n_dgn_sig_rm
eqtlgen_sig <- p_all %>% filter(p <= p_thre & module %in% module_use$module)
eqtlgen_sig_uniq <- eqtlgen_sig %>% distinct(meta) %>% pull()

eqtlgen_old_signal_uniq <- eqtlgen_old_signal %>% distinct(meta) %>% pull()


dim(dgn_sig)

sum(eqtlgen_all_snp_uniq %in% dgn_all_snp_uniq)

sum(dgn_sig_uniq %in% eqtlgen_sig_uniq)
sum(dgn_sig_uniq %in% eqtlgen_old_signal_uniq)
sum(dgn_sig_uniq %in% eqtlgen_all_snp_uniq)








dgn_sig %>%
  filter(module %in% paste0("module", module_use$module) ) %>%
  distinct(SNP_id)

dgn_sig %>%
  filter(SNP_id %in% eqtlgen_all_snp$SNP_id ) %>%
  distinct(SNP_id)


dgn_sig %>%
  filter(module %in% paste0("module", module_use$module) &
           SNP_id %in% eqtlgen_all_snp$SNP_id ) %>%
  distinct(SNP_id)




dgn_sig %>%
  filter(SNP_id %in% eqtlgen_sig_uniq) %>%
  left_join(eqtlgen_sig, by = c("SNP_id" = "meta") ) %>%
  group_by(SNP_id, module.x) %>%
  summarise("module_eqtlgen" = paste(sort(module.y), collapse = ";") ) %>%
  mutate("if_module_overlap" = module.x %in% paste0("module", str_split(module_eqtlgen, ";")[[1]]) ) %>%
  View()

  