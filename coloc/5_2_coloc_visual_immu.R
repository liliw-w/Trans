########## visualization of colocalization events ##########
########## by R package locuscomparer ##########
rm(list = ls())
library(data.table)
library(tidyverse)
library(locuscomparer)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

########## files and parameters ##########
gwas_pmid = as.numeric(args[1])
gwas_label = args[2]

pp4Thre = 0.75
popColoc = "EUR"
refGenom = 'hg19'

## folder for each gwas trait
dir_gwas = file.path("/scratch/midway2/liliw1/coloc",
                     paste0("pmid", gwas_pmid, "_", gwas_label))
dir_gwas_data = file.path(dir_gwas, "data")
dirVisData = file.path(dir_gwas, "colocVisData")
dirPlot = file.path(dir_gwas, "plot")

dir.create(dirVisData, showWarnings = FALSE)
dir.create(dirPlot, showWarnings = FALSE)

file_qtlColocReg_gwas = file.path(dir_gwas_data, "qtlColocReg_gwas.txt.gz")
file_gwasColocReg = file.path(dir_gwas_data, "gwasColocReg.txt.gz")
file_resColoc = file.path(dir_gwas_data, "resColoc.txt.gz")


########## read files ##########
qtlColocReg_gwas = fread(file_qtlColocReg_gwas)
gwasColocReg = fread(file_gwasColocReg)
resColoc = fread(file_resColoc)


########## only retain regions with high pp4 ##########
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

  file_qtlVis = file.path(dirVisData, paste0(reg, ".qtl.reg.coloc.txt") )
  file_gwasVis = file.path(dirVisData, paste0(reg, ".gwas.reg.coloc.txt") )

  fwrite(tmpqtlColocReg_gwas %>% select(rsid, Pval) %>% rename(pval = Pval), file_qtlVis, quote = FALSE, sep = "\t")
  fwrite(tmpgwasColocReg %>% select(rsid, pval), file_gwasVis, quote = FALSE, sep = "\t")

  figTitle = paste0("Region: ", tmpresColoc$Region, "; #SNPs: ", tmpresColoc$nsnps, "; PP.H4.abf: ", format(tmpresColoc$PP.H4.abf, digits = 4))
  fig <- locuscompare(in_fn1 = file_gwasVis, in_fn2 = file_qtlVis,
                      title1 = 'GWAS', title2 = 'QTL',
                      population = popColoc, genome = refGenom)

  fig <- fig + labs(title = figTitle) + theme(plot.title = element_text(hjust = 0.5, size = 12))

  ggsave(paste0(reg, ".reg.coloc.png"), fig, path = dirPlot, width = 6, height = 6)

  ### in progress
  cat(grep(reg, resColoc$Region, fixed = TRUE),
      "-th region:", reg, "(out of", nRegion, "regions)", "is done!", "\n")
}


########## Summary of the colocalized regions ##########
########## e.g. SNPs in each region, #regions for each module, #regions for each chr ##########
resColoc = resColoc %>% separate(col = Region, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)
resColoc$Module = factor(resColoc$Module, levels = paste0("module", 1:200))

## SNPs in each region
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

ggsave(paste0("coloc.region.summary.png"), fig, path = dirPlot, width = 10, height = 10)
