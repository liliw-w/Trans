##############################################
########### 50 msig pathways info ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(pathwayPCA)


# I/O & paras -----
file_msigdb_module <- '/project2/xuanyao/llw/MODULES/MSigDB/h.all.v7.4.symbols.gmt'

## output -----
file_qtl_table <- '/project2/xuanyao/llw/MODULES/MSigDB/msig_50_table.txt'


# read files -----
msigdb_module <- read_gmt(file_msigdb_module, description = TRUE)


# change col names and order -----
msigdb_module$size <- sapply(msigdb_module$pathways, length)

msigdb_module[c('TERMS', 'size', 'description', 'pathways')] %>%
  as_tibble() %>%
  arrange(desc(size)) %>%
  rename(
    pathway_genes = pathways
  ) %>%
  fwrite(
    file_qtl_table,
    quote = FALSE, sep = "\t", sep2 = c(""," ","")
  )

