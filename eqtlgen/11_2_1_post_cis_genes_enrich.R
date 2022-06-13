##############################################
########### enrichment of cis- genes ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
file_signal_cis_genes <- "postanalysis/signal_rm_infl_ratio_50.txt_cis_genes.txt"


# read files -----
signal_cis_genes <- fread(file_signal_cis_genes, header = TRUE)


# re-format to do enrichment -----
## g:profiler enrich -----
signal_cis_genes %>% distinct(nearest_gene) %>% pull() %>% paste(collapse = " ")

## DAVID enrich -----
signal_cis_genes %>% distinct(nearest_gene_symbol) %>% pull() %>% paste(collapse = " ")

