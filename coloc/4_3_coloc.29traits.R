rm(list = ls())
library(data.table)
library(tidyverse)
library(coloc)


file_pheno_manifest = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"
dir_gwas_data = "/project2/xuanyao/llw/GWAS/UKB_nealelab"
dir_coloc_data = "/project2/xuanyao/llw/coloc"
dir_coloc = "/scratch/midway2/liliw1/coloc"

pheno_manifest = fread(file_pheno_manifest)

########## files and parameters ##########
pop = "EUR"
gwasPhenocode_seq = pheno_manifest$phenocode_uniq
gwas_trait_type_seq = pheno_manifest$trait_type
gwas_trait_seq = pheno_manifest$trait
file_gwas_bgz_seq = file.path(dir_gwas_data, pheno_manifest$filename)
n_gwas = nrow(pheno_manifest)

#gwasPhenocode_seq = pheno_manifest$phenocode_uniq
p_included_thre = 1e-5

file_qtlColocReg = file.path(dir_coloc_data, "qtlColocReg.txt.gz")
file_snp_meta = file.path(dir_gwas_data, "full_variant_qc_metrics.txt.bgz")

snpMetaCol = c("chrom", "pos", "rsid", "high_quality")
snpMeta = fread(cmd = paste("gunzip -c", file_snp_meta), select = snpMetaCol)
qtlColocReg = fread(file_qtlColocReg)

for(i in 1:n_gwas){
  gwasPhenocode = gwasPhenocode_seq[i]
  gwas_trait_type = gwas_trait_type_seq[i]
  gwas_trait = gwas_trait_seq[i]
  file_gwas_bgz = file_gwas_bgz_seq[i]
  
  dir_gwas = file.path(dir_coloc,
                       paste("ukbb",
                             gwas_trait_type,
                             gwasPhenocode,
                             gwas_trait,
                             sep = "_")
  )
  dir.create(file.path(dir_gwas, "data"), recursive = TRUE, showWarnings = FALSE)
  
  
  ## output files
  setwd(dir_gwas)
  file_qtlColocReg_gwas = "data/qtlColocReg_gwas.txt.gz"
  file_gwasColocReg = "data/gwasColocReg.txt.gz"
  file_gwasRegTruncPthre = "data/gwasRegTruncPthre.txt"
  
  ########## read files ##########
  gwasCol = c("chr", "pos", paste(c("af", "beta", "se", "pval", "low_confidence"), pop, sep = "_") )
  gwas = fread(cmd = paste("gunzip -c", file_gwas_bgz), select = gwasCol)
  
  
  if(identical(gwas$pos, snpMeta$pos)){
    gwas = cbind(gwas, snpMeta[, c("rsid", "high_quality")])
  }else stop("GWAS variants and SNP meta are not aligned!")
  
  gwas$if_good = !is.na(gwas[[paste0("pval_", pop)]]) & gwas$high_quality & !gwas[[paste0("low_confidence_", pop)]]
  
  
  ########## gwas good SNPs & if overlapped with qtl SNPs & remove duplicated SNPs ##########
  gwas$if_qtl = paste(gwas$chr, gwas$pos, sep = ":") %in% qtlColocReg$SNP_ID
  gwas = gwas %>% filter(if_good & if_qtl &
                           !duplicated(paste(chr, pos, sep = ":")) )
  gwas$SNP_ID = paste(gwas$chr, gwas$pos, sep = ":")
  
  if(nrow(gwas) == 0) stop("No individuals in the specified population are phenotyped for this GWAS trait!")
  
  
  ########## qtl regions with overlapped SNPs ##########
  qtlColocReg_gwas = qtlColocReg[qtlColocReg$SNP_ID %in% gwas$SNP_ID, ]
  
  
  ########## coloc regions based on qtl ##########
  gwas$chr = as.integer(gwas$chr)
  
  gwasColocReg = qtlColocReg_gwas %>% select(c(Signal, Module, SNP_ID, Chr, Pos, Region)) %>% left_join(y = gwas, by = c("SNP_ID", "Chr"="chr", "Pos"="pos"))
  
  
  # range of pvalues in a region for GWAS and qtl & remove the regions whose min GWAS p < 1e-5
  gwasRegP = gwasColocReg %>% group_by(Region) %>% summarise(s = min(pval_EUR), l = max(pval_EUR))
  qtlRegP = qtlColocReg_gwas %>% group_by(Region) %>% summarise(s = min(Pval), l = max(Pval))
  
  gwasRegTruncPthre = gwasRegP$Region[gwasRegP$s < p_included_thre]
  
  
  ########## save results ##########
  fwrite(qtlColocReg_gwas, file_qtlColocReg_gwas, quote = FALSE, sep = "\t")
  fwrite(gwasColocReg, file_gwasColocReg, quote = FALSE, sep = "\t")
  fwrite(data.table(gwasRegTruncPthre), file_gwasRegTruncPthre, quote = FALSE, sep = "\t", col.names = FALSE)
  
  
  cat(i, '-th out of', n_gwas, "GWASs is done. Trait is", gwas_trait, ".\n")
  
  
  ##########################################
  #################  COLOC  ################
  ##########################################
  qtlN = 913
  qtlType = "quant"
  
  
  gwasType = ifelse(gwas_trait_type %in% c("continuous", "biomarkers"),
                    "quant",
                    "cc") # a trait is quantitative, trait_type is "continuous" or "biomarkers"
  
  pop = "EUR"
  pp4Thre = 0.75
  csThre = 0.95
  
  # output
  file_resColoc = "data/resColoc.txt.gz"
  file_resColocSNP = "data/resColocSNP.txt.gz"
  
  
  ########## files and parameter, read files ##########
  gwasRegTruncPthre = data.table("Region" = gwasRegTruncPthre)
  gwasTraitInfoCol = c("trait_type", "trait", "phenocode_uniq", "pheno_sex",
                       paste(c("n_cases", "n_controls"), pop, sep = "_") )
  gwasTraitInfo = fread(file_pheno_manifest, select = gwasTraitInfoCol)
  gwasTraitInfo = gwasTraitInfo %>% rename("phenocode" = "phenocode_uniq")
  
  ########## GWAS trait info, sample size
  gwasTraitInfo_usedGWAS = gwasTraitInfo %>% filter(phenocode == gwasPhenocode)
  if(nrow(gwasTraitInfo_usedGWAS) != 1) stop("There are more than one trait with the same given phenocode. Check more. e.g. The trait is not quantitative and have more than two phenotype categories. Or, e.g. The phenocode corresponds to multiple traits with different detailed descriptions.")
  n_cases = gwasTraitInfo_usedGWAS[[paste0("n_cases_", pop)]]
  n_controls = gwasTraitInfo_usedGWAS[[paste0("n_controls_", pop)]]
  gwasN = sum(n_cases, n_controls, na.rm = TRUE)
  
  
  ### remove SNPs whose se(\beta) equals 0, as coloc needs to use 1/se(\beta)
  if(any(gwasColocReg[[paste0("se_", pop)]] == 0)){
    gwasColocReg = gwasColocReg[gwasColocReg[[paste0("se_", pop)]] != 0, ]
    qtlColocReg_gwas = qtlColocReg_gwas[qtlColocReg_gwas$SNP_ID %in% gwasColocReg$SNP_ID, ]
  }
  
  
  ########## prepare coloc files and run coloc ##########
  resColoc = NULL
  resColocSNP = NULL
  nRegion = length(gwasRegTruncPthre$Region)
  for(reg in gwasRegTruncPthre$Region){
    ### extract the region for qtl and gwas
    tmpqtlColocReg_gwas = qtlColocReg_gwas %>% filter(Region == reg)
    tmpgwasColocReg = gwasColocReg %>% filter(Region == reg)
    
    
    ### construct coloc dataset D1 & D2
    D1 = list("pvalues" = tmpqtlColocReg_gwas$Pval,
              "N" = qtlN,
              "MAF" = tmpqtlColocReg_gwas$MAF,
              "type" = qtlType,
              "snp" = tmpqtlColocReg_gwas$SNP_ID)
    D2 = list("beta" = tmpgwasColocReg[[paste0("beta_", pop)]],
              "varbeta" = (tmpgwasColocReg[[paste0("se_", pop)]])^2,
              "type" = gwasType,
              "s" = if(gwasType=="cc") n_cases/gwasN else NULL,
              "snp" = tmpgwasColocReg$SNP_ID,
              "MAF" = tmpgwasColocReg[[paste0("af_", pop)]],
              "N" = gwasN)
    #"pvalues" = tmpgwasColocReg[[paste0("pval_", pop)]]
    if(gwasType=="cc"){
      D2 = within(D2, rm(MAF))
    }else{
      D2 = within(D2, rm(s))
    }
    
    ### do coloc
    coloc_res = coloc.abf(D1, D2)
    
    ### follow-up: credible set
    if(coloc_res$summary["PP.H4.abf"] > pp4Thre){
      o <- order(coloc_res$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(coloc_res$results$SNP.PP.H4[o])
      w <- which(cs > csThre)[1]
      resColocSNP = rbind(resColocSNP, data.table("Region" = reg,
                                                  "SNP_ID" = as.character(coloc_res$results[o,][1:w,]$snp),
                                                  "SNP.PP.H4" = coloc_res$results[o,][1:w,]$SNP.PP.H4))
    }
    
    ### organize result
    resColoc = rbind(resColoc, data.table("Region" = reg, t(coloc_res$summary)) )
    
    ### in progress
    cat(grep(reg, gwasRegTruncPthre$Region, fixed = TRUE),
        "-th region:", reg, "(out of", nRegion, "regions)", "is done!", "\n")
  }
  
  
  ########## Add p-value and trait info ##########
  resColoc = qtlColocReg %>% select(c(Signal, Pval)) %>%
    right_join(y = resColoc, by = c("Signal" = "Region")) %>%
    rename("Region" = "Signal") %>%
    mutate("phenocode" = gwasPhenocode)
  resColoc = gwasTraitInfo %>% select(c("phenocode", "trait")) %>%
    right_join(y = resColoc, by = c("phenocode")) %>%
    rename("Phenocode" = "phenocode")
  
  #resColoc = trait_info %>% select(c("GWAS ID", "Trait Abbreviation")) %>%
  #  right_join(y = resColoc, by = c("GWAS ID" = "Phenocode")) %>%
  #  rename("Phenocode" = "GWAS ID", "trait" = "Trait Abbreviation")
  
  
  ########## save results##########
  resColoc = resColoc %>% arrange(desc(PP.H4.abf))
  fwrite(resColoc, file_resColoc, quote = FALSE, sep = "\t")
  if(!is.null(resColocSNP)) fwrite(resColocSNP, file_resColocSNP, quote = FALSE, sep = "\t")
  
  
  cat(i, '-th out of', n_gwas, "GWASs is done. Trait is", gwas_trait, ".\n")
  
}
