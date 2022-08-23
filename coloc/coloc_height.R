########## Run coloc on regions for qtl and gwas traits ##########
########## These traits are ukbb traits from nealelab ##########
rm(list = ls())
library(data.table)
library(tidyverse)
library(coloc)


########## files and parameters ##########
#gwasPhenocode = 30080
qtlN = 913
qtlType = "quant"
gwasType = "quant"


#gwasType = ifelse(gwas_trait_type %in% c("continuous", "biomarkers"),
#                  "quant",
#                  "cc") # a trait is quantitative, trait_type is "continuous" or "biomarkers"
#af_


pop = "EUR"
pp4Thre = 0.75
csThre = 0.95

file_qtlColocReg_gwas = "qtlColocReg_gwas.txt.gz"
file_gwasColocReg = "gwasColocReg.txt.gz"
file_gwasRegTruncPthre = "gwasRegTruncPthre.txt"
#file_pheno_manifest = "/scratch/midway2/liliw1/UKB_nealelab/phenotype_manifest.tsv.bgz"
#file_pheno_manifest = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"

# Add p-value and trait info
file_qtlColocReg = "/project2/xuanyao/llw/coloc/qtlColocReg.txt.gz"
#file_trait_info = "/scratch/midway2/liliw1/coloc/ukbb_blood_traits.csv"

# output
file_resColoc = "resColoc.txt.gz"
file_resColocSNP = "resColocSNP.txt.gz"


########## files and parameter, read files ##########
qtlColocReg_gwas = fread(file_qtlColocReg_gwas)
gwasColocReg = fread(file_gwasColocReg)
gwasRegTruncPthre = fread(file_gwasRegTruncPthre, header = FALSE, col.names = "Region")

#gwasTraitInfoCol = c("trait_type", "trait", "phenocode_uniq", "pheno_sex",
#                     paste(c("n_cases", "n_controls"), pop, sep = "_") )
#gwasTraitInfo = fread(file_pheno_manifest, select = gwasTraitInfoCol)
#gwasTraitInfo = gwasTraitInfo %>% rename("phenocode" = "phenocode_uniq")

########## GWAS trait info, sample size
#gwasTraitInfo_usedGWAS = gwasTraitInfo %>% filter(phenocode == gwasPhenocode)
#if(nrow(gwasTraitInfo_usedGWAS) != 1) stop("There are more than one trait with the same given phenocode. Check more. e.g. The trait is not quantitative and have more than two phenotype categories. Or, e.g. The phenocode corresponds to multiple traits with different detailed descriptions.")
#n_cases = gwasTraitInfo_usedGWAS[[paste0("n_cases_", pop)]]
#n_controls = gwasTraitInfo_usedGWAS[[paste0("n_controls_", pop)]]
#gwasN = sum(n_cases, n_controls, na.rm = TRUE)


### remove SNPs whose se(\beta) equals 0, as coloc needs to use 1/se(\beta)
#if(any(gwasColocReg[[paste0("se_", pop)]] == 0)){
#  gwasColocReg = gwasColocReg[gwasColocReg[[paste0("se_", pop)]] != 0, ]
#  qtlColocReg_gwas = qtlColocReg_gwas[qtlColocReg_gwas$SNP_ID %in% gwasColocReg$SNP_ID, ]
#}


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
  D2 = list("pvalues" = tmpgwasColocReg$P,
            "N" = tmpgwasColocReg$N,
            "MAF" = tmpgwasColocReg$Freq_Tested_Allele_in_HRS,
            "beta" = tmpgwasColocReg$BETA,
            "varbeta" = (tmpgwasColocReg$SE)^2,
            "type" = gwasType,
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


########## Add p-value and trait info ##########
qtlColocReg = fread(file_qtlColocReg, header = TRUE)
#trait_info = fread(file_trait_info, sep = ",", header = TRUE)

resColoc = qtlColocReg %>% select(c(Signal, Pval)) %>%
  right_join(y = resColoc, by = c("Signal" = "Region")) %>%
  rename("Region" = "Signal") %>%
  mutate("Phenocode" = "height_Meta-analysis_Wood_et_al+UKBiobank_2018",
         "trait" = "height")

#resColoc = trait_info %>% select(c("GWAS ID", "Trait Abbreviation")) %>%
#  right_join(y = resColoc, by = c("GWAS ID" = "Phenocode")) %>%
#  rename("Phenocode" = "GWAS ID", "trait" = "Trait Abbreviation")


########## save results##########
resColoc = resColoc %>% arrange(desc(PP.H4.abf))
fwrite(resColoc, file_resColoc, quote = FALSE, sep = "\t")
if(!is.null(resColocSNP)) fwrite(resColocSNP, file_resColocSNP, quote = FALSE, sep = "\t")


