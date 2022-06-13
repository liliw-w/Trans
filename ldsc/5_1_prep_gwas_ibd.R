###############################################################
########### Prepare GWAS sum stats for trait ibd ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
gwasPhenocode <- 28067908
gwas_label <- "ibd"


file_gwas <- list.files(path = "/project2/xuanyao/llw/GWAS",
                        pattern = paste0("^pmid", gwasPhenocode, "_.*", gwas_label, ".*.txt.gz$"),
                        full.names = TRUE)
file_gwas_trait_info <- "/project2/xuanyao/llw/GWAS/72_traits_GWAS.csv"

file_snp_meta <- "/project2/xuanyao/llw/GWAS/UKB_nealelab/snp.meta.gwas.txt.gz"

file_out <- paste0("gwas/", gwasPhenocode, "_", gwas_label, ".tsv.gz")


# columns to be read from gwas file
gwasCol <- c("SNP", "Effect", "StdErr", "P.value", "Allele1", "Allele2")
gwasCol_newName <- c("SNP", "beta", "se", "P", "A1", "A2")


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

gwas <- snp_meta %>%
  select(SNP_ID, rsid) %>%
  right_join(gwas, by = c("SNP_ID" = "SNP")) %>%
  filter(!is.na(rsid) & !is.na(beta) & !duplicated(rsid))

if(nrow(gwas) == 0) stop("No SNPs left for this GWAS trait after QC!")

dim(gwas)


# 2. Capitalize alleles
gwas <- gwas %>% mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2))


# 3. calculate z-score -----
gwas$Z <- gwas$beta/gwas$se


# 4. add sample size N -----
gwasTraitInfo_usedGWAS <- gwas_trait_info %>% filter(Label == gwas_label & PMID == gwasPhenocode)
if(nrow(gwasTraitInfo_usedGWAS) != 1) stop("There are more than one trait with the same given phenocode. Check more. e.g. The trait is not quantitative and have more than two phenotype categories. Or, e.g. The phenocode corresponds to multiple traits with different detailed descriptions.")
gwas$N <- gwasTraitInfo_usedGWAS$N


# 5. rename columns and write out -----
gwas %>%
  select(rsid, N, Z, P, A1, A2) %>%
  rename(SNP = rsid) %>%
  fwrite(
    file = file_out,
    quote = FALSE,
    sep = "\t"
  )
