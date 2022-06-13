###############################################################
########### Prepare GWAS sum stats for trait ibd ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
gwasPhenocode <- 29083406
gwas_label <- "Allergy"


file_gwas <- list.files(path = "/project2/xuanyao/llw/GWAS",
                        pattern = paste0("^pmid", gwasPhenocode, "_.*", gwas_label, ".*.gz$"),
                        full.names = TRUE)
file_gwas_trait_info <- "/project2/xuanyao/llw/GWAS/72_traits_GWAS.csv"

file_snp_meta <- "/project2/xuanyao/llw/GWAS/UKB_nealelab/snp.meta.gwas.txt.gz"

file_out <- paste0("gwas/", gwasPhenocode, "_", gwas_label, ".tsv.gz")


# columns to be read from gwas file
gwasCol <- c("RS_ID", "BETA", "SE", "PVALUE", "EFFECT_ALLELE", "OTHER_ALLELE", "N")
gwasCol_newName <- c("rsid", "beta", "se", "P", "A1", "A2", "N")


# read data -----
gwas <- fread(
  file = file_gwas,
  select = gwasCol[!is.na(gwasCol)],
  col.names = gwasCol_newName,
  header = TRUE
)

gwas_trait_info <- fread(file_gwas_trait_info)
gwas_trait_info[gwas_trait_info == ""] <- NA

snp_meta <- fread(file_snp_meta)


# 1. append rs id to gwas snps and remove bad snps -----
dim(gwas)

gwas <- gwas %>%
  filter(!is.na(rsid) & !is.na(beta) & !duplicated(rsid))

if(nrow(gwas) == 0) stop("No SNPs left for this GWAS trait after QC!")

dim(gwas)


# 2. Capitalize alleles
gwas <- gwas %>% mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2))


# 3. calculate z-score -----
gwas$Z <- gwas$beta/gwas$se


# 5. rename columns and write out -----
gwas %>%
  select(rsid, N, Z, P, A1, A2) %>%
  rename(SNP = rsid) %>%
  fwrite(
    file = file_out,
    quote = FALSE,
    sep = "\t"
  )
