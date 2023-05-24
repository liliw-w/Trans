##############################################
########### table of signals of msig pathways in eQTLGen ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_qtl <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/signal_cis_genes_w_annot_rm_infl_ratio_50.txt'
file_msig_annot <- '/project2/xuanyao/llw/MODULES/MSigDB/msig_50_table.txt'

## output -----
file_qtl_table <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/qtl_table.txt'


# read files -----
qtl <- fread(file_qtl)
msig_annot <- fread(file_msig_annot, sep = '\t', sep2 = ' ', )


# change col names and order -----
qtl <- arrange(qtl, module, SNPChr, SNPPos) %>%
  mutate(
    P_adjusted_bonferroni = p*11*9918
  ) %>%
  left_join(
    select(msig_annot, TERMS, description_info),
    by = c('annot_module' = 'TERMS')
  ) %>%
  select(
    signal, module, annot_module, description_info, SNP, rsid, p, P_adjusted_bonferroni, 
    nearest_gene, near_genes, 
    gwas_catalog_trait, 
    gwas_catalog_report_gene, gwas_catalog_mapped_gene
  ) %>%
  rename(
    gene_module = module,
    module_msig = annot_module,
    module_description_info = description_info,
    P = p
  )

# merge gwas_catalog info for signals -----
qtl_gwas_uniq <- apply(qtl, 1, function(x){
  rbind(
    strsplit(x["gwas_catalog_trait"], ";") |> unlist() |> unique() |> paste(collapse = ";"),
    strsplit(x["gwas_catalog_report_gene"], ";") |> unlist() |> unique() |> paste(collapse = ";"),
    strsplit(x["gwas_catalog_mapped_gene"], ";") |> unlist() |> unique() |> paste(collapse = ";")
  )
}) %>%
  t() %>% as_tibble()
qtl[, c('gwas_catalog_trait', 'gwas_catalog_report_gene', 'gwas_catalog_mapped_gene')] <- qtl_gwas_uniq


# print out key message or write out -----
fwrite(
  qtl,
  file_qtl_table,
  quote = FALSE, sep = "\t"
)

