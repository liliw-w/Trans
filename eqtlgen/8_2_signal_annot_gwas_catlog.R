##############################################
########### annotate the signal pairs ###########
########### by gwas catalog snp-trait/disease for signal snps ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
ratio <- 50

file_signal_cis_genes <- paste0('postanalysis/signal_cis_genes_rm_infl_ratio_', ratio, '.txt')
file_gwas_catalog_hg19 <- '/project2/xuanyao/llw/gwas_catalog/gwas_catalog_hg19.txt'

dis_cis <- 1e+5

## output -----
file_signal_annot <- paste0('postanalysis/signal_cis_genes_w_annot_rm_infl_ratio_', ratio, '.txt')


# read files -----
signal_cis_genes <- fread(file_signal_cis_genes, header = TRUE)
gwas_catalog_hg19 <- fread(file_gwas_catalog_hg19, header = TRUE)


# organize data -----
## add (module:snp) column
signal_cis_genes <- mutate(signal_cis_genes, "signal" = paste0("module", module, ":", meta))



# add gwas catalog SNPs' relevant diseases/traits as annotation -----
gwas_catalog_hg19 <- mutate(
  gwas_catalog_hg19,
  "SNP_ID_hg19" = paste(CHR_ID_hg19, CHR_POS_hg19, sep = ":")
)


signal_gwas_catalog <- sapply(1:nrow(signal_cis_genes), function(x){
  tmp_signal_cis_genes = signal_cis_genes[x, ]
  
  ind_gwas = gwas_catalog_hg19$SNP_ID_hg19 %in% tmp_signal_cis_genes$meta
  
  tmp_cis_gwas_catalog_hg19 <- filter(gwas_catalog_hg19, CHR_ID_hg19 == !!tmp_signal_cis_genes$SNPChr) %>%
    mutate("dis" = abs(CHR_POS_hg19 - !!tmp_signal_cis_genes$SNPPos)) %>%
    filter(dis < !!dis_cis/2) %>%
    arrange(dis) %>%
    select(SNPS, dis, `DISEASE/TRAIT`, `REPORTED GENE(S)`, MAPPED_GENE)
  
  c(
    sum(ind_gwas) > 0,
    paste(unique(gwas_catalog_hg19[ind_gwas, SNPS]), collapse = ";"),
    paste(gwas_catalog_hg19[ind_gwas, `DISEASE/TRAIT`], collapse = ";"),
    paste(gwas_catalog_hg19[ind_gwas, `REPORTED GENE(S)`], collapse = ";"),
    paste(gwas_catalog_hg19[ind_gwas, MAPPED_GENE], collapse = ";"),
    
    nrow(tmp_cis_gwas_catalog_hg19) > 0,
    paste(tmp_cis_gwas_catalog_hg19$SNPS, collapse = ";"),
    paste(tmp_cis_gwas_catalog_hg19$dis, collapse = ";"),
    paste(tmp_cis_gwas_catalog_hg19$`DISEASE/TRAIT`, collapse = ";"),
    paste(tmp_cis_gwas_catalog_hg19$`REPORTED GENE(S)`, collapse = ";"),
    paste(tmp_cis_gwas_catalog_hg19$`MAPPED_GENE`, collapse = ";")
  )
  
})
signal_gwas_catalog <- as.data.table(t(signal_gwas_catalog))
colnames(signal_gwas_catalog) <- c(
  "if_gwas_catalog", "rsid", "gwas_catalog_trait",
  "gwas_catalog_report_gene", "gwas_catalog_mapped_gene",
  "if_cis_gwas_catalog", "cis_rsid", "cis_dis_gwas_catalog", "cis_gwas_catalog_trait",
  "cis_gwas_catalog_report_gene", "cis_gwas_catalog_mapped_gene"
)


# add the annotation columns to signal pairs & reorder rows and columns -----
signal_cis_genes <- bind_cols(signal_cis_genes, signal_gwas_catalog) %>%
  mutate(
    if_gwas_catalog = as.logical(if_gwas_catalog),
    if_cis_gwas_catalog = as.logical(if_cis_gwas_catalog),
    SNP = NULL
  ) %>%
  arrange(
    desc(if_gwas_catalog), desc(if_cis_gwas_catalog),
    module, SNPChr, SNPPos
  ) %>%
  relocate(
    signal,
    if_gwas_catalog, gwas_catalog_trait,
    if_cis_gwas_catalog, cis_gwas_catalog_trait,
    .before = everything()
  ) %>%
  rename(SNP = meta)


# print out key message or write out -----
fwrite(
  signal_cis_genes,
  file_signal_annot,
  quote = FALSE, sep = "\t"
)


