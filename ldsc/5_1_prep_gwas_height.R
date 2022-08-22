###############################################################
########### Prepare GWAS sum stats for trait ibd ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
file_gwas <- "/project2/xuanyao/llw/GWAS/height_Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz"

## output -----
file_out <- paste0("gwas/height_Meta-analysis_Wood_et_al+UKBiobank_2018.tsv.gz")


# read data -----
gwas <- fread(
  file = file_gwas,
  header = TRUE
)



# 1. remove missing snps -----
gwas <- gwas %>%
  filter(!is.na(SNP) & !is.na(BETA) & !duplicated(SNP))

if(nrow(gwas) == 0) stop("No SNPs left for this GWAS trait after QC!")


# 3. calculate z-score -----
gwas$Z <- gwas$BETA/gwas$SE


# 5. rename columns and write out -----
gwas %>%
  select(SNP, N, Z, P, Tested_Allele, Other_Allele) %>%
  rename(A1 = Tested_Allele, A2 = Other_Allele) %>%
  fwrite(
    file = file_out,
    quote = FALSE,
    sep = "\t"
  )
