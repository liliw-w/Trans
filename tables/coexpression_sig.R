##############################################
########### signals of co-expression modules ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_qtl <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/signal_cis_genes.txt'
file_enrich_top2 <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich_top2.txt'

## output -----
file_qtl_table <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/qtl_table.txt'


# read files -----
qtl <- fread(file_qtl)
enrich_top2 <- fread(file_enrich_top2)


# change col names and order -----
qtl <- separate(qtl, module, c(NA, "module"), sep = "module", convert = TRUE) %>%
  arrange(module, chr, pos) %>%
  rename(
    gene_module = module,
    P = p,
    FDR = q
  ) %>%
  mutate(
    chr = NULL,
    pos = NULL
  )

# add the top2 enrichment in the trans target modules -----
qtl <- left_join(
  qtl, enrich_top2,
  by = "gene_module"
)


# print out key message or write out -----
fwrite(
  qtl,
  file_qtl_table,
  quote = FALSE, sep = "\t"
)

