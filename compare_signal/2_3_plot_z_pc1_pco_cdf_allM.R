##############################################
########### plot cdf of abs z-scores of genes in a module for a signal SNP ###########
########### of all modules ###########
########### detected by methods only by trans-PCO, PC1, or both ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)

source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
# use already saved data from previous cdf plot
file_z_list <- list.files("pc1", "^M\\d+_z.rds$", full.names = TRUE)


# combine z's of signals of all modules & extract the max abs z for each module -----
df_z <- rbindlist(
  lapply(file_z_list, function(file_z){
    cat("File", file_z, "is running... \n\n")
    
    # read files -----
    tmp_df_z <- readRDS(file_z)
    
    # organize data -----
    module <- str_extract(basename(file_z), "\\d+") %>% as.numeric()
    
    # z for each snp across genes in a module
    mutate(tmp_df_z, "module" = !!module)
  })
)

# signal (snp, module) counts for each types: both PC1 and PCO, or only for PCO or PC1
df_sig <- df_z %>% group_by(type) %>% distinct(snp, module) %>% ungroup() %>% count(type)


# plot cdf -----
base_plt <- ggplot(df_z, aes(x = abs(z), color = type, fill = type)) +
  stat_ecdf(geom = "step", size = 0.8) +
  stat_ecdf(alpha = 0.2, geom = "area",color = NA) +
  labs(title = paste0("All modules"),
       x = quote(~"|Z|"), y = quote(~"F(|Z|)"),
       color = "Signal")

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

ggsave(paste0("pc1/M_all_z_cdf.pdf"), plt, height = 4, width = 4.5)
