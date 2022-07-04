##############################################
########### plot cdf of z-scores of signals detected by univariate and multivariate ###########
##############################################
rm(list = ls())
library(tidyverse)

source('~/Trans/plot/theme_my_pub.R')

# paras and I/O -----
file_df_z <- "uni/M_all_z.rds"


# read files -----
df_z <- readRDS(file_df_z)


# organize data for plotting -----
# number of signal snps for categories
df_sig <- count(df_z, type_trait)

# replace z inf values
inf_replace_val <- df_z %>% filter(!is.infinite(abs_z)) %>% pull(abs_z) %>% max() %>% ceiling()
df_z <- mutate(df_z, abs_z = replace(abs_z, is.infinite(abs_z), inf_replace_val))


# plot cdf -----
base_plt <- ggplot(df_z, aes(x = abs_z, linetype = type_trait)) +
  stat_ecdf(geom = "step", size = 1) +
  labs(x = quote(~"|Z|"), y = quote(~"F(|Z|)"), linetype = "Trait")

plt <- base_plt +
  scale_linetype_manual(
    breaks = c("Module", "Gene"),
    values = c("Module" = "solid", "Gene" = "dashed"),
    labels = c("Module" = paste0("Module (", filter(df_sig, type_trait == "Module") %>% pull(n), ")"),
               "Gene" = paste0("Gene (", filter(df_sig, type_trait == "Gene") %>% pull(n), ")")),
    guide = guide_legend(override.aes = list(size = 0.5))
  ) +
  theme_my_pub() +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    legend.background = element_blank(),
    legend.position = "bottom",
    #legend.key.size= unit(0.2, "cm"),
    legend.text = element_text(size = 10),
    
    #axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(paste0("uni/maxz_cdf.pdf"), plt, height = 4, width = 4.5)

