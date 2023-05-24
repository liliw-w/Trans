##############################################
### update to using regions whose lead-SNP are module-QTL signals
### previously, I used a pre-specified p threshold, e.g. 1e-8
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(cowplot)


# I/O & paras -----
pp4Thre = 0.75
pvalThre = 'module_QTL_sig'
nsnpsThre = 5

file_module_QTL_signals <- '/project2/xuanyao/llw/MODULES/MSigDB/FDR/signals.chr.module.perm10.fdr10.txt'
file_coexp_module = "/project2/xuanyao/llw/MODULES/MSigDB/result/coexp.module.rds"
dir_coloc_gwas <- list.files(
  "/scratch/midway2/liliw1/coloc_MSigDB",
  pattern = paste0("^ukbb_continuous_\\d+"),
  full.names = TRUE
)
file_resColoc <- sapply(
  dir_coloc_gwas,
  list.files,
  pattern = "^resColoc.txt.gz$", full.names = TRUE, recursive = TRUE
)


## output -----
dir_data_enrich = "ukbb_all/data_enrich"
file_resEnrich = paste0("ukbb_all/data_enrich/all.traits.enrich.", pvalThre, ".txt")


# read files -----
module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
coexp_module = readRDS(file_coexp_module)$moduleLabels
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
