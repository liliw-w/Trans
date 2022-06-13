###############################################################
########### Check overlapps of trans with cis ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


### I/O
file_trans <- '/scratch/midway2/liliw1/eQTGen_est_Sigma/postanalysis/LD.prun.in.chr.module.txt'
file_cis_e <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz'
file_cis_s <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_sQTL_signif_variant_gene_pairs.txt.gz'

file_trans_cis_e <- 'postanalysis/trans_cis_eQTL.txt'
file_trans_cis_s <- 'postanalysis/trans_cis_sQTL.txt'
file_sig_e <- 'postanalysis/signals_eQTL.txt'
file_sig_s <- 'postanalysis/signals_sQTL.txt'
file_gene_e <- 'postanalysis/Genes_eQTL.txt'
file_gene_s <- 'postanalysis/Genes_sQTL.txt'


### read sets
trans <- fread(file_trans, header = FALSE)
cis_e <- fread(file_cis_e, header = TRUE)
cis_s <- fread(file_cis_s, header = TRUE)


### overlaps between trans and cis
cis_e <- cis_e %>% filter(sid %in% trans$V1)
cis_s <- cis_s %>% filter(sid %in% trans$V1)


### output
# 1. trans_cis output
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


# 2. overlapped SNPs
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


# 3. corresponding e/s-Genes of overlapped SNPs
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


### message
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
  "cis-sQTL-Genes. \n"
)

