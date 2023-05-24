##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_qtl_rep <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/dgn_in_eqtlgen.txt'

## output -----
file_qtl_table  <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/dgn_in_eqtlgen_table.txt'


# read files -----
qtl_rep <- fread(file_qtl_rep)


# change col names and order -----
separate(qtl_rep, module, c(NA, "module"), sep = "module", convert = TRUE) %>%
  arrange(chr, pos, module) %>%
  rename(
    gene_module = module,
    SNP = SNP_id,
    P_DGN = p.dgn,
    P_eQTLGen = p.eqtlgen
  ) %>%
  mutate(
    chr = NULL,
    pos = NULL,
    q = NULL,
    if_rep = NULL
  ) %>%
  fwrite(
    file_qtl_table,
    quote = FALSE, sep = "\t"
  )

