##############################################
########### 166 Co-expression gene modules ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_coexp_module <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds"

## output -----
file_module_df <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/module.txt'


# read files -----
coexp_module <- readRDS(file_coexp_module)$moduleLabels


# organize data -----
df_module <- enframe(coexp_module, "gene", "gene_module") %>%
  group_by(gene_module) %>%
  summarise(
    module_size = n(),
    "genes_in_module" = paste(gene, collapse = ";")
  )

# print out key message or write out -----
fwrite(
  df_module,
  file_module_df,
  quote = FALSE, sep = "\t"
)

