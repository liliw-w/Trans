##############################################
########### z-scores across genes in a module for a signal SNP ###########
########### detected by methods only by trans-PCO, PC1, or both ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)

source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
dir_z <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/z/"
file_pco_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_pc1_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability_PC1/FDR/signals.chr.module.perm10.txt'


# read files -----
pco_sig <- fread(file_pco_sig, header = FALSE, col.names = c('signal', 'p', 'q'))
pc1_sig <- fread(file_pc1_sig, header = FALSE, col.names = c('signal', 'p', 'q'))


# organize data -----
rownames(pco_sig) <- pco_sig$signal
pco_sig <- pco_sig %>%
  separate(signal, into = c("module", "chr", "pos"), sep = "[:]", convert = TRUE) %>%
  separate(module, into = c(NA, "module"), sep = "module", convert = TRUE) %>%
  unite("snp", chr, pos, sep = ":", remove = FALSE)

rownames(pc1_sig) <- pc1_sig$signal
pc1_sig <- pc1_sig %>%
  separate(signal, into = c("module", "chr", "pos"), sep = "[:]", convert = TRUE) %>%
  separate(module, into = c(NA, "module"), sep = "module", convert = TRUE) %>%
  unite("snp", chr, pos, sep = ":", remove = FALSE)



# iterate across modules with signals by PC1 & extract z -----
sig_all <- union(select(pco_sig, snp, chr, module), select(pc1_sig, snp, chr, module))
module_seq <- pc1_sig %>% distinct(module) %>% pull()

for(module in module_seq){
  cat("Module", module, "is running... \n\n")
  
  # find chr's with signals by PC1
  chr_seq <- sig_all %>% filter(module == !!module) %>% distinct(chr) %>% pull()
  
  # combine signals across all chr's
  file_z <- paste0(dir_z, "z.module", module, ".chr", chr_seq, ".txt.gz")
  df_z <- lapply(file_z, function(x){
    cat(basename(x), "out of chr's:", chr_seq, "is running... \n\n")
    
    fread(x, header = TRUE) %>%
      filter(snp %in% !!sig_all$snp) %>%
      pivot_longer(-snp, names_to = "gene", values_to = "z")
  }) %>%
    rbindlist()
  
  
  # add signal type: if it's signal for both PC1 and PCO, or only for PCO or PC1 -----
  df_z <- df_z %>%
    mutate(
      "if_pco" = snp %in% (pco_sig %>% filter(module == !!module) %>% pull(snp)),
      "if_pc1" = snp %in% (pc1_sig %>% filter(module == !!module) %>% pull(snp)),
      "type" = paste(if_pco, if_pc1, sep = "_")
    ) %>%
    filter(type != "FALSE_FALSE")
  
  df_z$type <- factor(
    df_z$type,
    levels = c("TRUE_FALSE", "TRUE_TRUE", "FALSE_TRUE"),
    labels = c("Trans-PCO", "Both", "PC1")
  )
  
  
  # save z -----
  saveRDS(df_z, paste0("pc1/M", module, "_z.rds"))
}

