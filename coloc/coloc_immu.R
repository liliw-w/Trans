########## Run coloc on regions for qtl and gwas traits ##########
########## These traits are part of the 72 traits ##########
########## Focus on the 14 autoimmune diseases ##########
rm(list = ls())
library(data.table)
library(tidyverse)
library(coloc)

args = commandArgs(trailingOnly=TRUE)

########## files and parameters ##########
gwas_pmid_seq = c(24390342, 26192919, 26192919, 26192919, 26502338, 28067908, 28067908, 28067908, 29083406, 29892013, 30929738, 31604244)
gwas_label_seq = c("RA_GWASmeta_European", "CD", "IBD", "UC", "sle", "cd", "ibd", "uc", "Allergy", "AE", "ASTHMA", "MS")

qtlN = 913
qtlType = "quant"
gwasType = "cc"
pp4Thre = 0.75
csThre = 0.95

file_qtlColocReg = "/scratch/midway2/liliw1/coloc_MSigDB/qtlColocReg.txt.gz"
file_gwas_trait_info = "/project2/xuanyao/llw/GWAS/72_traits_GWAS.csv"

qtlColocReg = fread(file_qtlColocReg, header = TRUE)
gwas_trait_info = fread(file_gwas_trait_info)
gwas_trait_info[gwas_trait_info == ""] = NA

for(k in 1:length(gwas_pmid_seq)){
  gwas_pmid = gwas_pmid_seq[k]
  gwas_label = gwas_label_seq[k]
  
  ## folder for each gwas trait
  dir_gwas = file.path("/scratch/midway2/liliw1/coloc_MSigDB",
                       paste0("pmid", gwas_pmid, "_", gwas_label))
  dir_gwas_data = file.path(dir_gwas, "data")
  
  ## region files
  file_qtlColocReg_gwas = file.path(dir_gwas_data, "qtlColocReg_gwas.txt.gz")
  file_gwasColocReg = file.path(dir_gwas_data, "gwasColocReg.txt.gz")
  file_gwasRegTruncPthre = file.path(dir_gwas_data, "gwasRegTruncPthre.txt")
  
  ## files for adding p-value of region lead SNP and trait info
  
  ## output
  file_resColoc = file.path(dir_gwas_data, "resColoc.txt.gz")
  file_resColocSNP = file.path(dir_gwas_data, "resColocSNP.txt.gz")
 
  
  ########## read files ##########
  qtlColocReg_gwas = fread(file_qtlColocReg_gwas)
  gwasColocReg = fread(file_gwasColocReg)
  gwasRegTruncPthre = fread(file_gwasRegTruncPthre, header = FALSE, col.names = "Region")
  
  
  ########## GWAS trait info, sample size
  gwas_trait_info_usedGWAS = gwas_trait_info %>% filter(Label == gwas_label & PMID == gwas_pmid)
  n_cases = as.numeric(gwas_trait_info_usedGWAS[["N_case"]])
  n_controls = as.numeric(gwas_trait_info_usedGWAS[["N_control"]])
  gwasN = sum(n_cases, n_controls)
  
  ## remove SNPs whose se(\beta) equals 0, as coloc needs to use 1/se(\beta)
  if(any(gwasColocReg[["se"]] == 0)){
    gwasColocReg = gwasColocReg[gwasColocReg[["se"]] != 0, ]
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
    D2 = list("pvalues" = tmpgwasColocReg[["pval"]],
              "N" = if( !is.na(gwasN) ) gwasN else tmpgwasColocReg[["N"]],
              "MAF" = tmpgwasColocReg[["af"]],
              "beta" = tmpgwasColocReg[["beta"]],
              "varbeta" = (tmpgwasColocReg[["se"]])^2,
              "type" = gwasType,
              "s" = if(gwasType=="cc" & !is.na(gwasN) ) n_cases/gwasN else if(gwasType=="cc" & is.na(gwasN) ) tmpgwasColocReg[["s"]] else NULL,
              "snp" = tmpgwasColocReg$SNP_ID)
    
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
  
  
  ########## Add p-value and trait info to coloc regions ##########
  resColoc = qtlColocReg %>% select(c(Signal, Pval)) %>%
    right_join(y = resColoc, by = c("Signal" = "Region")) %>%
    rename("Region" = "Signal") %>%
    mutate("gwas_label" = gwas_label)
  resColoc = gwas_trait_info %>%
    filter(PMID == gwas_pmid & Label == gwas_label) %>%
    select(c("Trait", "Label")) %>%
    right_join(y = resColoc, by = c("Label" = "gwas_label")) %>%
    rename("trait" = "Trait", "Phenocode" = "Label")
  
  
  ########## save results##########
  resColoc = resColoc %>% arrange(desc(PP.H4.abf))
  fwrite(resColoc, file_resColoc, quote = FALSE, sep = "\t")
  if(!is.null(resColocSNP)) fwrite(resColocSNP, file_resColocSNP, quote = FALSE, sep = "\t")
  
  str(resColoc)
  str(sum(resColoc$PP.H4.abf > pp4Thre))
  
  cat("Trait", gwas_label, ". \n")
  
}








