##############################################
########### plot pdf of z-scores of categorized signals detected by univariate and multivariate ###########
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

# number of signal snps for categories
df_sig <- df_z %>%
  group_by(type_trait, type) %>%
  summarise("abs_z" = max(abs_z), "n" = n()) %>%
  ungroup() %>%
  mutate("type_trait_num" = as.numeric(type_trait))


# for pairwise p comparison
#my_comparisons <- list( c("Gene", "Module") )

# plot violin & boxplot & stat_compare_means -----
base_plt <- ggplot(df_z, aes(x = type, y = abs_z)) +
  geom_violin(
    aes(color = type, fill = type,
        linetype = type_trait, size = type_trait),
    position = position_dodge(width = 0.5)
  ) +
  stat_summary(
    aes(group = type_trait),
    fun = mean,
    geom = "point",
    shape = 18,
    size = 2,
    color = "black",
    show.legend = FALSE,
    position = position_dodge(width = 0.5)
  ) +
  #stat_compare_means(comparisons = my_comparisons,
  #                   #label.y = y_max + c(-0.2, 0, 0.5),
  #                   vjust = 2, show.legend = FALSE) +
  geom_text(
    data = df_sig,
    aes(label = n, color = type, group = type_trait),
    position = position_dodge(width = 0.5),
    vjust = -0.1, size = 3
  ) +
  labs(x = NULL, y = quote(~"|Z|"), linetype = "Trait", color = "Signal")

plt <- base_plt +
  scale_linetype_manual(
    breaks = c("Module", "Gene"),
    values = c("Module" = "solid", "Gene" = "twodash"),
    guide = guide_legend(override.aes = list(size = 0.3))
  ) +
  scale_size_manual(
    breaks = c("Module", "Gene"),
    values = c("Module" = 0.5, "Gene" = 0.3),
    guide = "none"
  ) +
  scale_colour_manual(
    breaks = c("Multi", "Both", "Uni"),
    values = c("Multi" = "#6a1424", "Both" = "#002080", "Uni" = "#e89c31"),
    guide = "none"
  ) +
  scale_fill_manual(
    breaks = c("Multi", "Both", "Uni"),
    values = c("Multi" = "#c28c96", "Both" = "#ccccff", "Uni" = "#f8e1c1"),
    guide = "none"
  ) +
  theme_my_pub() +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    legend.background = element_blank(),
    legend.position = "bottom",
    legend.key.size= unit(0.2, "cm"),
    legend.text = element_text(size = 10)
  )
plt

ggsave(paste0("uni/maxz_cat_sig_pdf.pdf"), plt, height = 3.5, width = 4.5)

