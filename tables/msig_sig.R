##############################################
########### table of signals of coexpression modules in eQTLGen ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_qtl <- '/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/signal_cis_genes_w_annot.txt'
file_msig_annot <- '/project2/xuanyao/llw/MODULES/MSigDB/msig_50_table.txt'

## output -----
file_qtl_table <- '/project2/xuanyao/llw/MODULES/MSigDB/postanalysis/qtl_table.txt'


# read files -----
qtl <- fread(file_qtl)
msig_annot <- fread(file_msig_annot, sep = '\t', sep2 = ' ', )


# change col names and order -----
arrange(qtl, module, chr, pos) %>%
  left_join(
    select(msig_annot, TERMS, description_info),
    by = c('annot_module' = 'TERMS')
  ) %>%
  select(
    signal, module, annot_module, description_info, 
    SNP, rsid, p, q, 
    nearest_gene, near_genes, 
    #gwas_catalog_trait, 
    #gwas_catalog_report_gene, gwas_catalog_mapped_gene
  ) %>%
  rename(
    gene_module = module,
    module_msig = annot_module,
    module_description_info = description_info,
    P = p,
    FDR = q
  ) %>%
  fwrite(
    file_qtl_table,
    quote = FALSE, sep = "\t"
  )

