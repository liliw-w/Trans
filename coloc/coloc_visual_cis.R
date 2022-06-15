########## visualization of colocalization events ##########
########## by R package locuscomparer ##########
rm(list = ls())
library(data.table)
library(tidyverse)
library(locuscomparer)
library(cowplot)


## output files
dir_gwas = '/scratch/midway2/liliw1/coloc_cis_ARHGEF3/'
setwd(dir_gwas)

file_qtlColocReg_gwas = "qtlColocReg_gwas.txt.gz"
file_gwasColocReg = "gwasColocReg.txt.gz"
file_resColoc_all = "resColoc.txt.gz"

pp4Thre = 0.75
popColoc = "EUR"
refGenom = 'hg19'


########## files and parameter, read files ##########
qtlColocReg_gwas = fread(file_qtlColocReg_gwas)
gwasColocReg = fread(file_gwasColocReg)
resColoc_all = fread(file_resColoc_all)


########## only retain regions with high pp4 ##########
#removeReg = c("module4:12:8142264")
#resColoc = resColoc %>% filter(PP.H4.abf > pp4Thre)


## remove SNPs without rsid, as locuscomparer requires rsid
qtlColocReg_gwas = qtlColocReg_gwas[qtlColocReg_gwas$rsid != "", ]
gwasColocReg$rsid <- qtlColocReg_gwas[match(gwasColocReg$SNP_ID, qtlColocReg_gwas$SNP_ID), ]$rsid
gwasColocReg = gwasColocReg[gwasColocReg$rsid != "", ]


gwasPhenocode_seq = unique(resColoc_all$trait)
gwas_trait_seq = unique(resColoc_all$trait)
n_gwas = length(unique(resColoc_all$trait))


# i <- 57, 196, 264, 272, 396:398 ??? "There must be one and only one chromosome."
for(i in 1:n_gwas){
  gwasPhenocode = gwasPhenocode_seq[i]
  gwas_trait = gwas_trait_seq[i]
  
  dir.create(file.path("colocVisData"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("plot"), recursive = TRUE, showWarnings = FALSE)


  resColoc <- resColoc_all %>% filter(trait == gwasPhenocode)
  
  ########## Visualize each region & save ##########
  nRegion = length(resColoc$Region)
  
  for(reg in resColoc$Region){
    tmpresColoc = resColoc %>% filter(Region == reg)
    tmpqtlColocReg_gwas = qtlColocReg_gwas %>% filter(Region == reg)
    tmpgwasColocReg = gwasColocReg %>%
      filter(gene == gwasPhenocode &
               Region == reg &
               SNP_ID %in% tmpqtlColocReg_gwas$SNP_ID)
    
    file_qtlVis = paste0("colocVisData/", reg, ".", gwasPhenocode, ".qtl.reg.coloc.txt")
    file_gwasVis = paste0("colocVisData/", reg, ".", gwasPhenocode, ".gwas.reg.coloc.txt")
    
    fwrite(tmpqtlColocReg_gwas %>% select(rsid, Pval) %>% rename(pval = Pval),
           file_qtlVis, quote = FALSE, sep = "\t")
    fwrite(tmpgwasColocReg %>% select(rsid, npval) %>% rename(pval = npval),
           file_gwasVis, quote = FALSE, sep = "\t")
    
    figTitle = paste0("Region: ", tmpresColoc$Region,
                      "; Gene: ", gwasPhenocode,
                      "; #SNPs: ", tmpresColoc$nsnps,
                      "; PP.H4.abf: ", format(tmpresColoc$PP.H4.abf, digits = 4))
    fig <- locuscompare(in_fn1 = file_gwasVis, in_fn2 = file_qtlVis,
                        title1 = 'Gene', title2 = 'QTL',
                        population = popColoc, genome = refGenom)
    
    fig <- fig + labs(title = figTitle) + theme(plot.title = element_text(hjust = 0.5, size = 12))
    
    ggsave(paste0(reg, ".", gwasPhenocode, ".reg.coloc.png"), fig, path = "plot/", width = 6, height = 6)
    
    ### in progress
    cat(grep(reg, resColoc$Region, fixed = TRUE),
        "-th region:", reg, "(out of", nRegion, "regions)", "is done!", "\n")
    cat("Gene: ", gwasPhenocode, '\n')
  }
}


### Look at how many regions have coloc cis- genes; coloc to how many cis- genes; what are these genes
nsnpsThre <- 5
pvalThre <- "module_QTL"

res_coloc_reg_prop = resColoc_all %>%
  group_by(Region) %>%
  summarise(nGene = n(),
            nGeneColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) ) %>%
  ungroup()
res_coloc_reg_prop <- resColoc_all %>%
  filter(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) %>%
  group_by(Region) %>%
  summarise(gene = paste(trait, collapse = ";") ) %>%
  right_join(res_coloc_reg_prop, by = "Region") %>%
  arrange(desc(nGeneColoc), desc(nGene) )

res_coloc_reg_prop = res_coloc_reg_prop %>%
  mutate(propColoc = nGeneColoc/nGene )


fwrite(res_coloc_reg_prop,  paste0("coloc_region_prop_pvalThre-", pvalThre, ".txt"),
       quote = FALSE, sep = "\t")






