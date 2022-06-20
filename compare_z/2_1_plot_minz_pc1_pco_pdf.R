##############################################
########### plot pdf (violin plot) of min abs z-scores across genes in a module for a signal SNP ###########
########### detected by methods only by trans-PCO, PC1, or both ###########
##############################################
rm(list = ls())
library(tidyverse)
library(ggpubr)

source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
# use already saved data from previous cdf plot
file_z_list <- list.files("pc1", "^M\\d+_z.rds$", full.names = TRUE)


for(file_z in file_z_list){
  cat("File", file_z, "is running... \n\n")
  
  # read files -----
  df_z <- readRDS(file_z)
  
  # organize data -----
  module <- str_extract(basename(file_z), "\\d+")
  
  # min abs z for each snp across genes in a module
  df_z <- df_z %>%
    pivot_wider(names_from = gene, values_from = z) %>%
    rowwise(snp:type) %>%
    summarise("min_abs_z" = min(abs(c_across(where(is.numeric))))) %>%
    ungroup()
  df_sig <- df_z %>% count(type)
  
  # for pairwise p comparison
  my_comparisons <- list( c("Trans-PCO", "Both"), c("Both", "PC1"), c("Trans-PCO", "PC1") )
  y_max <- max(abs(df_z$min_abs_z))
  
  
  # plot violin & boxplot & stat_compare_means -----
  base_plt <- ggplot(df_z, aes(x = type, y = min_abs_z, color = type, fill = type)) +
    geom_violin() +
    #geom_quasirandom(varwidth = FALSE, shape = 16, size = 1, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black") +
    stat_compare_means(comparisons = my_comparisons,
                       vjust = 2, show.legend = FALSE) +
    #stat_compare_means() +
    labs(title = paste0("M", module),
         x = NULL, y = quote(~"min|Z|"))
  
  plt <- base_plt +
    scale_x_discrete(
      breaks = c("Trans-PCO", "Both", "PC1"),
      labels = c("Trans-PCO" = paste0("Trans-PCO \n (", filter(df_sig, type == "Trans-PCO") %>% pull(n), ")"),
                 "Both" = paste0("Both \n (", filter(df_sig, type == "Both") %>% pull(n), ")"),
                 "PC1" = paste0("PC1 \n (", filter(df_sig, type == "PC1") %>% pull(n), ")") )
    ) +
    scale_colour_manual(
      breaks = c("Trans-PCO", "Both", "PC1"),
      values = c("Trans-PCO" = "#6a1424", "Both" = "#002080", "PC1" = "#d08c2c"),
      guide = "none"
    ) +
    scale_fill_manual(
      breaks = c("Trans-PCO", "Both", "PC1"),
      values = c("Trans-PCO" = "#c28c96", "Both" = "#ccccff", "PC1" = "#f8e1c1"),
      guide = "none"
    ) +
    theme_my_pub() +
    theme(
      #panel.grid.major.y = element_line(linetype = "dashed"),
      
      legend.background = element_blank(),
      legend.position = "bottom",
      legend.key.size= unit(0.2, "cm"),
      legend.text = element_text(size = 10),
      
      plot.title = element_text(vjust = -1),
      
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  
  ggsave(paste0("pc1/M", module, "_minz_pdf.pdf"), plt, height = 4, width = 5)
  
}
