##############################################
########### extract gwas catalog SNPs into a bed file to do liftover ###########
########### from hg38 build to hg19 build ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_gwas_catalog <- '/project2/xuanyao/llw/gwas_catalog/gwas_catalog_v1.0.2-associations_e107_r2022-08-26.tsv'

## output -----
file_gwas_catalog_snp_hg38 <- "/project2/xuanyao/llw/gwas_catalog/gwas_catalog_snp_hg38.bed"


# read files -----
gwas_catalog <- fread(file_gwas_catalog, header = TRUE, quote = "")



# extract unique SNPs, remove SNPs in haplotype or interaction, add word chr for a bed format -----
gwas_catalog_snp_hg38 <- select(
  gwas_catalog,
  CHR_ID, CHR_POS
) %>%
  distinct(
    CHR_ID, CHR_POS, .keep_all = TRUE
  ) %>%
  mutate(
    CHR_ID = as.numeric(CHR_ID),
    CHR_POS = as.numeric(CHR_POS)
  ) %>%
  filter(complete.cases(.)) %>%
  mutate(
    "CHR_ID" = paste0("chr", CHR_ID),
    "end" = CHR_POS
  )


# print out key message or write out -----
fwrite(
  gwas_catalog_snp_hg38,
  file_gwas_catalog_snp_hg38,
  col.names = FALSE, quote = FALSE, sep = "\t"
)

