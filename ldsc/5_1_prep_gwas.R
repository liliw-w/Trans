###############################################################
########### Prepare GWAS sum stats ###########
###############################################################
#rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
gwasPhenocode <- as.numeric(snakemake@params[['gwasPhenocode']])

dir_gwas_data <- "/project2/xuanyao/llw/GWAS/UKB_nealelab"
pop <- "EUR"

file_gwas_bgz <- list.files(dir_gwas_data, paste0(".*-", gwasPhenocode, "-.*.tsv.bgz"), full.names = TRUE)
file_snp_meta <- list.files(dir_gwas_data, "snp.meta.gwas.txt.gz", full.names = TRUE)
file_gwasTraitInfo <- list.files(dir_gwas_data, "phenotype_manifest.tsv.bgz", full.names = TRUE)

file_out <- paste0("gwas/", gwasPhenocode, ".tsv.gz")

#file_pheno_manifest <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv"
#pheno_manifest <- fread(file_pheno_manifest)
#gwasPhenocode = pheno_manifest$phenocode_uniq[1]
#gwas_trait_type = pheno_manifest$trait_type[1]
#gwas_trait = pheno_manifest$trait[1]

# read data -----
gwasCol <- c("chr", "pos", "ref", "alt", paste(c("pval", "beta", "se", "low_confidence"), pop, sep = "_") )
gwas <- fread(cmd = paste("gunzip -c", file_gwas_bgz), select = gwasCol)

snp_meta <- fread(file_snp_meta)

gwasTraitInfoCol <- c("trait_type", "phenocode", "pheno_sex", "description",
                      paste(c("n_cases", "n_controls"), pop, sep = "_") )
gwasTraitInfo <- fread(cmd = paste("gunzip -c", file_gwasTraitInfo), select = gwasTraitInfoCol)


# 1. append rs id to gwas snps -----
if(identical(gwas$pos, snp_meta$pos)){
  gwas <- cbind(gwas, snp_meta[, c("rsid", "high_quality")])
}else stop("GWAS variants and SNP meta are not aligned!")


# 2. remove bad snps -----
gwas$if_good <- !is.na(gwas[[paste0("pval_", pop)]]) &
  gwas$high_quality &
  !gwas[[paste0("low_confidence_", pop)]]
gwas <- gwas %>% filter(if_good & !duplicated(paste(chr, pos, sep = ":")) )


# 3. calculate z-score -----
gwas$Z <- gwas$beta/gwas$se


# 4. add sample size N -----
gwasTraitInfo_usedGWAS <- gwasTraitInfo %>% filter(phenocode == gwasPhenocode)
if(nrow(gwasTraitInfo_usedGWAS) != 1) stop("There are more than one trait with the same given phenocode. Check more. e.g. The trait is not quantitative and have more than two phenotype categories. Or, e.g. The phenocode corresponds to multiple traits with different detailed descriptions.")
gwas$N <- sum(gwasTraitInfo_usedGWAS[[paste0("n_cases_", pop)]], gwasTraitInfo_usedGWAS[[paste0("n_controls_", pop)]], na.rm = TRUE)


# check if any errors -----
if(nrow(gwas) == 0) stop("No individuals in the specified population are phenotyped for this GWAS trait!")


# 5. rename columns and write out -----
gwas %>%
  select(rsid, N, Z, starts_with("pval"), alt, ref) %>%
  rename(A1 = alt,
         A2 = ref,
         SNP = rsid,
         P = starts_with("pval")) %>%
  fwrite(
    file = file_out,
    quote = FALSE,
    sep = "\t"
  )
