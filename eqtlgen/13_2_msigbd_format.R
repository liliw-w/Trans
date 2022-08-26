##############################################
########### re-format MSigDB gmt format to data frame ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)
library(pathwayPCA)


# I/O & paras -----
file_msigdb_module <- '/project2/xuanyao/llw/MODULES/MSigDB/h.all.v7.4.symbols.gmt'

## output -----
file_msigdb_module_df <- '/project2/xuanyao/llw/MODULES/MSigDB/h.all.v7.4.symbols.txt'


# read files -----
msigdb_module <- read_gmt(file_msigdb_module, description = TRUE)


# convert gmt format to data frame -----
names(msigdb_module$pathways) <- msigdb_module$TERMS
msigdb_module <- bind_rows(
  lapply(msigdb_module$pathways, function(x) tibble("gene" = x)),
  .id = "category"
)


# print out key message and write out -----
fwrite(
  msigdb_module,
  file_msigdb_module_df,
  quote = FALSE, sep = "\t"
)
