###############################################################
########### Check overlapps of trans with cis ###########
########### trans signals are from both coexpression modules & msig pathways ###########
########### cis signals are from DGN cis-e & cis-s ###########
########### but it makes more sense to use eQTLGen cis-e ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O -----
file_trans1 <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/signal_uniq_rm_infl_ratio_50.txt'
file_trans2 <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/signal_uniq_rm_infl_ratio_50.txt'
file_cis_e <- '/project2/xuanyao/llw/DGN_data/DGN_sumstat/DGN_eQTL/DGN_eQTL_signif_variant_gene_pairs.txt.gz'
file_cis_s <- '/project2/xuanyao/llw/DGN_data/DGN_sumstat/DGN_sQTL/DGN_sQTL_signif_variant_gene_pairs.txt.gz'

## output -----
file_sig_e <- 'trans_coexp_msig_cis_dgn/signals_eQTL.txt'
file_sig_s <- 'trans_coexp_msig_cis_dgn/signals_sQTL.txt'
file_trans_cis_e <- 'trans_coexp_msig_cis_dgn/trans_cis_eQTL.txt'
file_trans_cis_s <- 'trans_coexp_msig_cis_dgn/trans_cis_sQTL.txt'
file_gene_e <- 'trans_coexp_msig_cis_dgn/Genes_eQTL.txt'
file_gene_s <- 'trans_coexp_msig_cis_dgn/Genes_sQTL.txt'


# read sets -----
trans1 <- fread(file_trans1, header = TRUE)
trans2 <- fread(file_trans2, header = TRUE)
trans <- rbind(trans1, trans2) %>% select(meta) %>% distinct()
cis_e <- fread(file_cis_e, header = TRUE)
cis_s <- fread(file_cis_s, header = TRUE)


# overlaps between trans and cis -----
cis_e <- cis_e %>% filter(sid %in% !!trans$meta)
cis_s <- cis_s %>% filter(sid %in% !!trans$meta)


# output -----
## 1. overlapped SNPs
fwrite(
  cis_e %>% distinct(sid),
  file_sig_e,
  sep = "\t", quote = FALSE, col.names = FALSE
)
fwrite(
  cis_s %>% distinct(sid),
  file_sig_s,
  sep = "\t", quote = FALSE, col.names = FALSE
)


## 2. trans_cis output
fwrite(
  cis_e,
  file_trans_cis_e,
  sep = "\t", quote = FALSE
)
fwrite(
  cis_s,
  file_trans_cis_s,
  sep = "\t", quote = FALSE
)


## 3. corresponding e/s-Genes of overlapped SNPs
fwrite(
  cis_e %>% distinct(pid),
  file_gene_e,
  sep = "\t", quote = FALSE, col.names = FALSE
)
fwrite(
  cis_s %>% distinct(pid),
  file_gene_s,
  sep = "\t", quote = FALSE, col.names = FALSE
)


# message -----
cat(
  cis_e %>% distinct(sid) %>% nrow(),
  "trans-eQTLs (out of",
  trans %>% nrow(),
  ") are also cis-eQTLs, corresponding to",
  cis_e %>% distinct(pid) %>% nrow(),
  "cis-eQTL-Genes. \n",
  
  "\n \n",
  
  cis_s %>% distinct(sid) %>% nrow(),
  "trans-eQTLs (out of",
  trans %>% nrow(),
  ") are also cis-sQTLs, corresponding to",
  cis_s %>% distinct(pid) %>% nrow(),
  "cis-sQTL-Genes. \n",
  
  "\n \n",
  
  rbind(cis_e, cis_s) %>% distinct(sid) %>% nrow(),
  "trans-eQTLs (out of",
  trans %>% nrow(),
  ") are either cis-eQTLs or cis-sQTLs. \n"
)

