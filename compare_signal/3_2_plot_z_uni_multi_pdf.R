##############################################
########### plot pdf of z-scores of signals detected by univariate and multivariate ###########
##############################################
rm(list = ls())
library(tidyverse)
library(ggpubr)

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


# for pairwise p comparison
my_comparisons <- list( c("Gene", "Module") )


# plot violin & boxplot & stat_compare_means -----
base_plt <- ggplot(df_z, aes(x = type_trait, y = abs_z)) +
  geom_violin(aes(linetype = type_trait)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     #label.y = y_max + c(-0.2, 0, 0.5),
                     vjust = 2, show.legend = FALSE) +
  labs(x = NULL, y = quote(~"|Z|"), linetype = "Trait")


plt <- base_plt +
  scale_linetype_manual(
    breaks = c("Module", "Gene"),
    values = c("Module" = "solid", "Gene" = "twodash"),
    labels = c("Module" = paste0("Module (", filter(df_sig, type_trait == "Module") %>% pull(n), ")"),
               "Gene" = paste0("Gene (", filter(df_sig, type_trait == "Gene") %>% pull(n), ")")),
    guide = guide_legend(override.aes = list(size = 0.3))
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

ggsave(paste0("uni/maxz_pdf.pdf"), plt, height = 4, width = 4.5)

