##############################################
########### covert gwas catalog snps to hg19 build ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_gwas_catalog <- '/project2/xuanyao/llw/gwas_catalog/gwas_catalog_v1.0.2-associations_e107_r2022-08-26.tsv'
file_gwas_catalog_snp_hg19 <- '/project2/xuanyao/llw/gwas_catalog/gwas_catalog_snp_hg19.bed'


## output -----
file_gwas_catalog_hg19 <- '/project2/xuanyao/llw/gwas_catalog/gwas_catalog_hg19.txt'


# read files -----
gwas_catalog <- fread(file_gwas_catalog, header = TRUE, quote = "")
gwas_catalog_snp_hg19 <- fread(
  file_gwas_catalog_snp_hg19,
  header = FALSE,
  select = c(1, 2, 4),
  col.names = c("CHR_ID_hg19", "CHR_POS_hg19", "snp_hg38")
)


# organize data -----
## reformat gwas catalog snps of hg19 to make it match with original file -----
gwas_catalog_snp_hg19 <- separate(
  gwas_catalog_snp_hg19,
  snp_hg38,
  c("CHR_ID", NA, "CHR_POS"),
  sep = "[:-]"
) %>%
  separate(
    CHR_ID_hg19,
    c(NA, "CHR_ID_hg19"),
    sep = "chr"
  ) %>%
  separate(
    CHR_ID,
    c(NA, "CHR_ID"),
    sep = "chr"
  )


# add hg19 meta to gwas catalog file -----
gwas_catalog <- left_join(
  gwas_catalog, gwas_catalog_snp_hg19,
  by = c("CHR_ID", "CHR_POS")
)


# update SNPs' rsid with SNP's current rsid -----
gwas_catalog[MERGED == 1, "SNPS"] <- paste0("rs", gwas_catalog[MERGED == 1, SNP_ID_CURRENT])


# print out key message or write out -----
fwrite(
  select(
    gwas_catalog,
    `DISEASE/TRAIT`,
    SNPS, CHR_ID_hg19, CHR_POS_hg19, CHR_ID, CHR_POS,
    `REPORTED GENE(S)`, `MAPPED_GENE`
  ),
  file_gwas_catalog_hg19,
  quote = FALSE, sep = "\t"
)
