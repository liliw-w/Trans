########## visualization of colocalization events ##########
########## by R package locuscomparer ##########
rm(list = ls())
library(data.table)
library(tidyverse)
library(locuscomparer)
library(cowplot)

file_pheno_manifest = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"
dir_coloc = "/scratch/midway2/liliw1/coloc"

pheno_manifest = fread(file_pheno_manifest)

########## files and parameters ##########
gwasPhenocode_seq = pheno_manifest$phenocode_uniq
gwas_trait_type_seq = pheno_manifest$trait_type
gwas_trait_seq = pheno_manifest$trait
n_gwas = nrow(pheno_manifest)


########## files and parameter ##########
pp4Thre = 0.75
popColoc = "EUR"
refGenom = 'hg19'

for(i in 1:n_gwas){
  gwasPhenocode = gwasPhenocode_seq[i]
  gwas_trait_type = gwas_trait_type_seq[i]
  gwas_trait = gwas_trait_seq[i]
  
  dir_gwas = file.path(dir_coloc,
                       paste("ukbb",
                             gwas_trait_type,
                             gwasPhenocode,
                             gwas_trait,
                             sep = "_")
  )
  dir.create(file.path(dir_gwas, "colocVisData"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_gwas, "plot"), recursive = TRUE, showWarnings = FALSE)
  
  
  ## output files
  setwd(dir_gwas)
  file_qtlColocReg_gwas = "data/qtlColocReg_gwas.txt.gz"
  file_gwasColocReg = "data/gwasColocReg.txt.gz"
  file_resColoc = "data/resColoc.txt.gz"
  
  ########## files and parameter, read files ##########
  qtlColocReg_gwas = fread(file_qtlColocReg_gwas)
  gwasColocReg = fread(file_gwasColocReg)
  resColoc = fread(file_resColoc)
  
  
  ########## only retain regions with high pp4 ##########
  #removeReg = c("module4:12:8142264")
  resColoc = resColoc %>% filter(PP.H4.abf > pp4Thre)
  
  
  ########## Visualize each region & save ##########
  nRegion = length(resColoc$Region)
  
  ## remove SNPs without rsid, as locuscomparer requires rsid
  qtlColocReg_gwas = qtlColocReg_gwas[qtlColocReg_gwas$rsid != "", ]
  gwasColocReg = gwasColocReg[gwasColocReg$rsid != "", ]
  
  for(reg in resColoc$Region){
    tmpresColoc = resColoc %>% filter(Region == reg)
    tmpqtlColocReg_gwas = qtlColocReg_gwas %>% filter(Region == reg)
    tmpgwasColocReg = gwasColocReg %>% filter(Region == reg)
    
    file_qtlVis = paste0("colocVisData/", reg, ".qtl.reg.coloc.txt")
    file_gwasVis = paste0("colocVisData/", reg, ".gwas.reg.coloc.txt")
    
    fwrite(tmpqtlColocReg_gwas %>% select(rsid, Pval) %>% rename(pval = Pval),
           file_qtlVis, quote = FALSE, sep = "\t")
    fwrite(tmpgwasColocReg %>% select(rsid, starts_with("pval")) %>% rename(pval = pval_EUR),
           file_gwasVis, quote = FALSE, sep = "\t")
    
    figTitle = paste0("Region: ", tmpresColoc$Region, "; #SNPs: ", tmpresColoc$nsnps, "; PP.H4.abf: ", format(tmpresColoc$PP.H4.abf, digits = 4))
    fig <- locuscompare(in_fn1 = file_gwasVis, in_fn2 = file_qtlVis,
                        title1 = 'GWAS', title2 = 'QTL',
                        population = popColoc, genome = refGenom)
    
    fig <- fig + labs(title = figTitle) + theme(plot.title = element_text(hjust = 0.5, size = 12))
    
    ggsave(paste0(reg, ".reg.coloc.png"), fig, path = "plot/", width = 6, height = 6)
    
    ### in progress
    cat(grep(reg, resColoc$Region, fixed = TRUE),
        "-th region:", reg, "(out of", nRegion, "regions)", "is done!", "\n")
  }
  
  
  ########## Summary of the colocalized regions ##########
  ########## e.g. SNPs in each region, #regions for each module, #regions for each chr ##########
  resColoc = resColoc %>% separate(col = Region, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)
  resColoc$Module = factor(resColoc$Module, levels = paste0("module", 1:200))
  
  ##SNPs in each region
  figNumofSnps <- ggplot(resColoc, aes(y = reorder(Region, nsnps), x = nsnps)) +
    geom_bar(stat = "identity") +
    labs(x = "#SNPs in region", y = "Region") +
    theme_bw() + theme(axis.text.y = element_text(size = 7), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # module frequency
  figModCount <- ggplot(resColoc, aes(y = Module)) +
    geom_bar(stat = "count") +
    labs(x = "#regions") +
    theme_bw() + theme(axis.text.y = element_text(size = 7), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # chr frequency
  figChrCount <- ggplot(resColoc, aes(y = Chr)) +
    geom_bar(stat = "count") +
    labs(x = "#regions") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # combine the plots
  right_col <- plot_grid(figModCount, figChrCount, labels = c('B', "C"), ncol = 1)
  fig <- plot_grid(figNumofSnps, right_col, labels = c('A', ''), ncol = 2)
  
  ggsave("coloc.region.summary.png", fig, path = "plot/", width = 10, height = 10)
  
  
  cat(i, '-th out of', n_gwas, "GWASs is done. Trait is", gwas_trait, ".\n")
  
}

