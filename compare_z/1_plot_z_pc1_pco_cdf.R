##############################################
########### plot cdf of z-scores across genes in a module for a signal SNP ###########
########### detected by methods only by trans-PCO, PC1, or both ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggpubr)

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


sig_all <- union(select(pco_sig, snp, chr, module), select(pc1_sig, snp, chr, module))

module_seq <- pc1_sig %>% distinct(module) %>% pull()


for(module in module_seq){
  cat("Module", module, "is running... \n\n")
  
  chr_seq <- sig_all %>% filter(module == !!module) %>% distinct(chr) %>% pull()
  
  
  file_z <- paste0(dir_z, "z.module", module, ".chr", chr_seq, ".txt.gz")
  
  df_z <- lapply(file_z, function(x){
    cat(basename(x), "out of chr's:", chr_seq, "is done. \n")
    
    fread(x, header = TRUE) %>%
      filter(snp %in% !!sig_all$snp) %>%
      pivot_longer(-snp, names_to = "gene", values_to = "z")
  }) %>%
    rbindlist()
  
  
  # add signal type: if it's signal for both PC1 and PCO, or only for PCO or PC1 -----
  df_z <- df_z %>%
    mutate(
      "if_pco" = snp %in% pco_sig$snp,
      "if_pc1" = snp %in% pc1_sig$snp,
      "type" = paste(if_pco, if_pc1, sep = "_")
    )
  df_z$type <- factor(
    df_z$type,
    levels = c("TRUE_FALSE", "TRUE_TRUE", "FALSE_TRUE"),
    labels = c("Trans-PCO", "Both", "PC1")
  )
  
  
  # number of signal snps for categories of both PC1 and PCO, or only for PCO or PC1
  df_sig <- df_z %>% group_by(type) %>% distinct(snp) %>% ungroup() %>% count(type)
  
  
  # 3. cdf
  base_plt <- ggplot(df_z, aes(x = abs(z), color = type, fill = type)) +
    stat_ecdf(geom = "step", size = 1) +
    stat_ecdf(alpha = 0.2, geom = "area",color = NA) +
    labs(title = paste0("M", module),
         x = quote(~"|Z|"), y = quote(~"F(|Z|)"),
         color = "Method")
  
  plt <- base_plt +
    scale_colour_manual(
      breaks = c("Trans-PCO", "Both", "PC1"),
      values = c("Trans-PCO" = "#6a1424", "Both" = "#002080", "PC1" = "#e89c31"),
      labels = c("Trans-PCO" = paste0("Trans-PCO (", filter(df_sig, type == "Trans-PCO") %>% pull(n), ")"),
                 "Both" = paste0("Both (", filter(df_sig, type == "Both") %>% pull(n), ")"),
                 "PC1" = paste0("PC1 (", filter(df_sig, type == "PC1") %>% pull(n), ")") ),
      guide = guide_legend(override.aes = list(size = 1))
    ) +
    scale_fill_manual(
      breaks = c("Trans-PCO", "Both", "PC1"),
      values = c("Trans-PCO" = "#c28c96", "Both" = "#ccccff", "PC1" = "#eeb96e"),
      guide = "none"
    ) +
    theme_my_pub() +
    theme(
      panel.grid.major.y = element_line(linetype = "dashed"),
      
      legend.background = element_blank(),
      legend.position = "bottom",
      legend.key.size= unit(0.2, "cm"),
      legend.text = element_text(size = 10),
      
      #axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      plot.title = element_text(vjust = -1),
      
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  
  saveRDS(plt, paste0("pc1/M", module, "_z.rds"))
  ggsave(paste0("pc1/M", module, "_z_cdf.pdf"), plt, height = 4.5, width = 4.5)
  
}


