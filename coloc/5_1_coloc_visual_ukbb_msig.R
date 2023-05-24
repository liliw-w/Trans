##############################################
########## visualization of colocalization events ##########
########## by R package locuscomparer ##########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(locuscomparer)
library(cowplot)


# loop over all blood traits -----
file_blood_trait_info <- "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv"
trait_info <- fread(file_blood_trait_info, sep = ",", header = TRUE)
pheno_seq <- trait_info$`GWAS ID`
n_gwas <- length(pheno_seq)

for(gwasPhenocode in pheno_seq){
  # I/O & paras -----
  pp4Thre <- 0.75
  popColoc <- "EUR"
  refGenom <- 'hg19'
  
  
  ## coloc directory for the trait
  dir_coloc_gwas <- list.files(
    "/project2/xuanyao/llw/MODULES/MSigDB/coloc_MSigDB",
    pattern = paste0("ukbb_continuous_", gwasPhenocode),
    full.names = TRUE
  )
  dir.create(file.path(dir_coloc_gwas, "colocVisData"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_coloc_gwas, "plot"), recursive = TRUE, showWarnings = FALSE)
  
  file_qtlColocReg_gwas <- paste(dir_coloc_gwas, "data/qtlColocReg_gwas.txt.gz", sep = '/')
  file_gwasColocReg <- paste(dir_coloc_gwas, "data/gwasColocReg.txt.gz", sep = '/')
  file_resColoc <- paste(dir_coloc_gwas, "data/resColoc.txt.gz", sep = '/')
  
  
  ## output -----
  file_plt_coloc_sum <- paste(dir_coloc_gwas, "plot/coloc.region.summary.pdf", sep = '/')
  
  
  # read files -----
  qtlColocReg_gwas <- fread(file_qtlColocReg_gwas)
  gwasColocReg <- fread(file_gwasColocReg)
  resColoc <- fread(file_resColoc)
  
  ## only retain regions with high pp4
  resColoc <- resColoc %>% filter(PP.H4.abf > pp4Thre)
  
  
  # 1. Visualize each region -----
  nRegion <- length(resColoc$Region)
  
  ## remove SNPs without rsid, as locuscomparer requires rsid
  qtlColocReg_gwas <- qtlColocReg_gwas[qtlColocReg_gwas$rsid != "", ]
  gwasColocReg <- gwasColocReg[gwasColocReg$rsid != "", ]
  
  for(reg in resColoc$Region){
    file_qtlVis <- paste0(dir_coloc_gwas, "/colocVisData/", reg, ".qtl.reg.coloc.txt")
    file_gwasVis <- paste0(dir_coloc_gwas, "/colocVisData/", reg, ".gwas.reg.coloc.txt")
    file_plt_vis <- paste0(dir_coloc_gwas, "/plot/", reg, ".reg.coloc.pdf")
    
    
    ## input files of locuscomparer -----
    tmpresColoc <- resColoc %>% filter(Region == reg)
    tmpqtlColocReg_gwas <- qtlColocReg_gwas %>% filter(Region == reg)
    tmpgwasColocReg <- gwasColocReg %>% filter(Region == reg)
    
    fwrite(tmpqtlColocReg_gwas %>% select(rsid, Pval) %>% rename(pval = Pval),
           file_qtlVis, quote = FALSE, sep = "\t")
    fwrite(tmpgwasColocReg %>% select(rsid, starts_with("pval")) %>% rename(pval = pval_EUR),
           file_gwasVis, quote = FALSE, sep = "\t")
    
    
    ## locus vis -----
    figTitle <- paste0("Region: ", tmpresColoc$Region, "; #SNPs: ", tmpresColoc$nsnps, "; PP.H4.abf: ", format(tmpresColoc$PP.H4.abf, digits = 4))
    fig <- locuscompare(in_fn1 = file_gwasVis, in_fn2 = file_qtlVis,
                        title1 = 'GWAS', title2 = 'QTL',
                        population = popColoc, genome = refGenom)
    fig <- fig + labs(title = figTitle) + theme(plot.title = element_text(hjust = 0.5, size = 12))
    
    
    ## save vis -----
    ggsave(file_plt_vis, fig, width = 4, height = 2.5)
    
    cat(grep(reg, resColoc$Region, fixed = TRUE),
        "-th region:", reg, "(out of", nRegion, "regions)", "is done!", "\n")
  }
  
  
  # 2. Summary of the colocalized regions -----
  resColoc <- resColoc %>% separate(col = Region, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)
  resColoc$Module <- factor(resColoc$Module, levels = paste0("module", 1:200))
  
  ## SNPs in each region -----
  figNumofSnps <- ggplot(resColoc, aes(y = reorder(Region, nsnps), x = nsnps)) +
    geom_bar(stat = "identity") +
    labs(x = "#SNPs in region", y = "Region") +
    theme_bw() + theme(axis.text.y = element_text(size = 7), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ## number of regions for each module -----
  figModCount <- ggplot(resColoc, aes(y = Module)) +
    geom_bar(stat = "count") +
    labs(x = "#regions") +
    theme_bw() + theme(axis.text.y = element_text(size = 7), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ## number of regions for each chr -----
  figChrCount <- ggplot(resColoc, aes(y = Chr)) +
    geom_bar(stat = "count") +
    labs(x = "#regions") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ## combine the plots -----
  right_col <- plot_grid(figModCount, figChrCount, labels = c('B', "C"), ncol = 1)
  fig <- plot_grid(figNumofSnps, right_col, labels = c('A', ''), ncol = 2)
  
  
  ## save vis -----
  ggsave(file_plt_coloc_sum, fig, width = 4, height = 5)
  
  cat(which(gwasPhenocode == pheno_seq), '-th out of', n_gwas, "GWASs is done. \n\n")
  
}

