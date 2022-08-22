rm(list = ls())
library(data.table)
library(tidyverse)
library(cowplot)

setwd('/project2/xuanyao/llw/coloc/immune_traits/')

########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 'module_QTL_sig'
nsnpsThre = 5
#gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)

file_resColoc = list.files(".", "resColoc.txt.gz", recursive = TRUE)

dirPlot = "~/xuanyao_llw/coloc/immune_traits/pmid_all"
file_plot = file.path(dirPlot, paste0("all.traits.coloc.region.summary.pvalThre-", pvalThre, ".png"))
file_res_coloc_reg_prop = file.path(dirPlot, paste0("coloc_region_prop_pvalThre-", pvalThre, ".txt"))

gwas_pmid_seq = c(29892013, 31604244, 24390342, 29083406, 30929738, 26502338, 26192919, 26192919, 26192919, 28067908, 28067908, 28067908)
gwas_label_seq = c("AE", "MS", "RA_GWASmeta_European", "Allergy", "ASTHMA", "sle", "IBD", "CD", "UC", "ibd", "cd", "uc")


file_module_QTL_signals = '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
range(module_QTL_signals$p)
range(module_QTL_signals$q)


########## loop over all traits ##########
########## number of regions & coloc regions, proportion of coloc regions ##########
resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

resColoc_all = resColoc_all %>% mutate('if_module_QTL' = Region %in% module_QTL_signals$signal )

res_coloc_reg_prop = resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegion = n(),
            nRegionColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre),
            nRegionPval = sum(if_module_QTL),
            nRegionPvalColoc = sum(if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)) %>%
  ungroup()
res_coloc_reg_prop = res_coloc_reg_prop %>%
  mutate(propColoc = nRegionColoc/nRegion, propPvalColoc = nRegionPvalColoc/nRegionPval) %>%
  arrange(desc(propPvalColoc), desc(propColoc))


########## save results ##########
fwrite(res_coloc_reg_prop, file_res_coloc_reg_prop, quote = FALSE, sep = "\t")
