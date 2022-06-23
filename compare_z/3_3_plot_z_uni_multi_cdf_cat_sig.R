##############################################
########### plot cdf of z-scores of categorized signals detected by univariate and multivariate ###########
##############################################
rm(list = ls())
library(tidyverse)
library(ggrepel)


source('~/Trans/plot/theme_my_pub.R')

# paras and I/O -----
file_df_z <- "uni/M_all_z.rds"


# read files -----
df_z <- readRDS(file_df_z)


# organize data for plotting -----
# replace z inf values
inf_replace_val <- df_z %>% filter(!is.infinite(abs_z)) %>% pull(abs_z) %>% max() %>% ceiling()
df_z <- mutate(df_z, abs_z = replace(abs_z, is.infinite(abs_z), inf_replace_val))

# number of signal snps for categories & annotation labels
df_sig <- df_z %>%
  group_by(type_trait, type) %>%
  summarise("abs_z" = max(abs_z), "y" = 1, "n" = n()) %>%
  ungroup() %>%
  mutate("type_trait_num" = as.numeric(type_trait))


# plot cdf -----
base_plt <- ggplot(df_z, aes(x = abs_z, linetype = type_trait, color = type)) +
  stat_ecdf(geom = "step", size = 1) +
  geom_text_repel(data = df_sig, aes(y = y, label = n, segment.linetype = type_trait_num),
                  nudge_y = -0.1, direction = "x", hjust = "right", size = 3) +
  labs(x = quote(~"|Z|"), y = quote(~"F(|Z|)"), linetype = "Trait", color = "Signal")

plt <- base_plt +
  scale_linetype_manual(
    breaks = c("Module", "Gene"),
    values = c("Module" = "solid", "Gene" = "twodash"),
    guide = guide_legend(override.aes = list(size = 0.5))
  ) +
  scale_colour_manual(
    breaks = c("Multi", "Both", "Uni"),
    values = c("Multi" = "#6a1424", "Both" = "#002080", "Uni" = "#e89c31"),
    guide = guide_legend(override.aes = list(label = ""))
  ) +
  theme_my_pub() +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    legend.background = element_blank(),
    legend.position = "right",
    #legend.key.size= unit(0.2, "cm"),
    legend.text = element_text(size = 10),
    
    #axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(paste0("uni/maxz_cat_sig_cdf.pdf"), plt, height = 3.5, width = 5)


