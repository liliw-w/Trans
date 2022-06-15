### update to using regions whose lead-SNP are module-QTL signals
### previously, I used a pre-specified p threshold, e.g. 1e-8

setwd('/project2/xuanyao/llw/coloc/immune_traits')

rm(list = ls())
library(data.table)
library(tidyverse)
library(cowplot)

file_module_QTL_signals = '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
range(module_QTL_signals$p)
range(module_QTL_signals$q)


########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 'module_QTL_sig'
nsnpsThre = 5
#gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)

file_coexp_module = "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds"
dir_data_enrich = "pmid_all/data_enrich"
file_resEnrich = paste0("pmid_all/data_enrich/all.traits.enrich.", pvalThre, ".txt")

coexp_module = readRDS(file_coexp_module)$moduleLabels


########## loop over all traits, region v.s. module ##########
#file_resColoc = paste0("/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/pheno", gwasPhenocode_seq, ".resColoc.txt.gz")
file_resColoc = list.files(".", "resColoc.txt.gz", recursive = TRUE)
resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

resColoc_all = resColoc_all %>% mutate('if_module_QTL' = Region %in% module_QTL_signals$signal )


resEnrich = resColoc_all %>%
  filter(if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) %>%
  select(c(Phenocode, trait, Region, Pval, nsnps, PP.H4.abf)) %>%
  separate(col = Region, into = c("Module", NA, NA), sep = ":", remove = FALSE, convert = TRUE)


########## Add genes for module ##########
coexp_module = data.table("Module" = paste0("module", coexp_module), "Gene" = names(coexp_module))
resEnrich$query = paste0(">", apply(resEnrich, 1, function(x) paste(x, collapse = ";") ))
resEnrich$gene = sapply(resEnrich$Module, function(x) paste(coexp_module[coexp_module$Module == x, Gene], collapse = " ") )


########## write out in format of g:profile multiple query ##########
res = resEnrich %>% select(c(Phenocode, trait, query, gene)) %>%
  pivot_longer(c(query, gene), names_to = NULL, values_to = "enrich")
res %>% group_by(Phenocode, trait) %>%
  group_walk(~ fwrite(.x,
                      file.path(dir_data_enrich,
                                paste0("pheno", .y$Phenocode, ".", .y$trait, ".enrich.", pvalThre, ".txt")),
                      quote = FALSE, col.names = FALSE) )

fwrite(resEnrich, file_resEnrich, quote = FALSE, sep = "\t")
