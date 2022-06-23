##############################################
########### z-scores of signals detected by univariate and multivariate ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)

# paras and I/O -----
dir_z <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/z/"
file_pco_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
file_uni_sig <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_transeQTL_battle2014.txt'

file_df_z <- "uni/M_all_z.rds"


# read files -----
pco_sig <- fread(file_pco_sig, header = FALSE, col.names = c('signal', 'p', 'q'))
uni_sig <- fread(file_uni_sig, header = TRUE)


# organize data -----
pco_sig <- pco_sig %>%
  separate(signal, into = c("module", "chr", "pos"), sep = "[:]", convert = TRUE) %>%
  separate(module, into = c(NA, "module"), sep = "module", convert = TRUE) %>%
  unite("snp", chr, pos, sep = ":", remove = FALSE)


# abs z of (SNP, gene) of uni signals -----
uni_sig <- uni_sig %>%
  mutate("PVAL" = 10^(-LOG_PVAL),
         'abs_z' = sqrt(qchisq(p = PVAL, df = 1, lower.tail = FALSE)))


# max abs z of (SNP, module) signal pairs of PCO signals across modules and chr -----
module_seq <- pco_sig %>% distinct(module) %>% pull()

maxabsz_pco <- lapply(module_seq, function(module){
  cat("Module", module, "is running... \n\n")
  
  # chr's with signals
  sig_module = pco_sig %>% filter(module == !!module)
  chr_seq = distinct(sig_module, chr) %>% pull()
  
  # combine signals across all chr's
  file_z = paste0(dir_z, "z.module", module, ".chr", chr_seq, ".txt.gz")
  lapply(file_z, function(x){
    cat(basename(x), "out of chr's:", chr_seq, "is running... \n\n")
    
    fread(x, header = TRUE) %>%
      filter(snp %in% !!sig_module$snp) %>%
      rowwise(snp) %>%
      summarise("max_abs_z" = max(abs(c_across(where(is.numeric))))) %>%
      ungroup() %>%
      mutate("module" = paste0("M", !!module))
  }) %>%
    rbindlist()
  
}) %>%
  rbindlist()


# rbind uni and multi signals & add signal type -----
# if it's signal for both uni and multi, or only for uni or multi
df_z <- bind_rows(
  uni_sig %>%
    select(GENE_NAME, snp, abs_z) %>%
    rename(trait = GENE_NAME) %>%
    mutate("if_multi" = snp %in% !!maxabsz_pco$snp,
           "if_uni" = TRUE,
           "type_trait" = "Gene"),
  
  maxabsz_pco %>%
    rename(trait = module, abs_z = max_abs_z) %>%
    mutate("if_multi" = TRUE,
           "if_uni" = snp %in% !!uni_sig$snp,
           "type_trait" = "Module")
) %>%
  mutate("type" = paste(if_multi, if_uni, sep = "_"))

df_z$type <- factor(
  df_z$type,
  levels = c("TRUE_FALSE", "TRUE_TRUE", "FALSE_TRUE"),
  labels = c("Multi", "Both", "Uni")
)

df_z$type_trait <- factor(
  df_z$type_trait,
  levels = c("Module", "Gene")
)


# save z -----
saveRDS(df_z, file_df_z)
