##############################################
########### ld clumping file prep for gwas ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)

# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

pop <- "EUR"

file_pheno_manifest <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv"
file_snp_meta <- "/project2/xuanyao/llw/GWAS/UKB_nealelab/full_variant_qc_metrics.txt.bgz"

## output -----
file_unique_snp_prefix <- 'ld_clump/uniq_snp_ukb_nealelab'
file_with_p_prefix <- 'ld_clump/ld_clump_snp_withp'


# read files -----
pheno_manifest <- data.table::fread(file_pheno_manifest)

snpMetaCol <- c("chrom", "pos", "rsid", "high_quality")
snpMeta <- data.table::fread(cmd = paste("gunzip -c", file_snp_meta), select = snpMetaCol)


# unique snps -----
snpMeta$SNP <- paste(snpMeta$chrom, snpMeta$pos, sep = ":")

filter(snpMeta, high_quality) %>%
  distinct(chrom, SNP) %>%
  group_by(chrom) %>%
  group_walk(~ data.table::fwrite(
    .x, 
    file = str_glue('{file_unique_snp_prefix}_chr{.y$chrom}.txt'), 
    sep = '\t', quote = FALSE, col.name = FALSE
  ))




# gwas trait snp with p -----
for(k in seq(nrow(pheno_manifest))){
  ## determine trait -----
  gwas_trait_type <- "continuous"
  gwasPhenocode <- pheno_manifest$`GWAS ID`[k]
  gwas_trait <- pheno_manifest$`Trait Abbreviation`[k]
  
  
  ## determine output -----
  file_with_p <- str_glue('{file_with_p_prefix}_{gwas_trait_type}-{gwasPhenocode}.txt.gz')
  
  
  ## read gwas files -----
  file_gwas_bgz <- list.files(
    "/project2/xuanyao/llw/GWAS/UKB_nealelab",
    paste0("^", gwas_trait_type, "-", gwasPhenocode, ".*.tsv.bgz$"),
    full.names = TRUE
  )
  if(length(file_gwas_bgz) == 0) stop("No GWAS trait!")
  
  gwasCol <- c("chr", "pos", paste(c("af", "beta", "se", "pval", "low_confidence"), pop, sep = "_") )
  gwasCol_new <- c("chr", "pos", "af", "beta", "se", "pval", "low_confidence")
  gwas <- data.table::fread(
    cmd = paste("gunzip -c", file_gwas_bgz), 
    select = gwasCol, 
    col.names = gwasCol_new
  )
  
  
  ## SNP filtering for gwas -----
  if(identical(gwas$pos, snpMeta$pos)){
    gwas <- cbind(gwas, snpMeta[, c("rsid", "high_quality")])
  }else stop("GWAS variants and SNP meta are not aligned!")
  
  gwas <- gwas[gwas$pval < 5e-8, ]
  
  gwas$if_good <- !is.na(gwas$pval) & gwas$high_quality & !gwas$low_confidence
  gwas$SNP_ID <- paste(gwas$chr, gwas$pos, sep = ":")
  gwas <- gwas %>% filter(if_good & !duplicated(SNP_ID))
  
  
  ## snp with p -----
  select(gwas, SNP_ID, pval) %>% 
    rename(SNP = SNP_ID, P = pval) %>% 
    data.table::fwrite(
      ., 
      file = file_with_p,
      sep = '\t', quote = FALSE
    )
  
  cat(str_glue('{k}-th trait is done. \n'))
}

