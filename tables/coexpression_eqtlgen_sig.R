##############################################
########### table of signals of coexpression modules in eQTLGen ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_qtl <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/signal_cis_genes_w_annot_rm_infl_ratio_50.txt'
file_enrich_top2 <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich_top2.txt'

## output -----
file_qtl_table <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/qtl_table.txt'


# read files -----
qtl <- fread(file_qtl)
enrich_top2 <- fread(file_enrich_top2)

# change fdr level to 0.5 to be consistent with numbers in manuscript
qtl <- filter(qtl, p < 0.05/129/9918)


# change col names and order -----
qtl <- arrange(qtl, module, SNPChr, SNPPos) %>%
  mutate(
    P_adjusted_bonferroni = p*129*9918
  ) %>%
  select(
    signal, module, SNP, rsid, p, P_adjusted_bonferroni, 
    nearest_gene, near_genes, 
    gwas_catalog_trait, 
    gwas_catalog_report_gene, gwas_catalog_mapped_gene
  ) %>%
  rename(
    gene_module = module,
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

