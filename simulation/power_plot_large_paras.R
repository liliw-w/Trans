###############################################################
########### Power plot when paras are large ###########
###############################################################
rm(list = ls())

library(tidyverse)
library(ggbeeswarm)

# files & paras -----
file_power_all <- '/project2/xuanyao/llw/simulation_lambda0.1/new_Sigma/power.large.paras.rds'
xlab_name <- "Sample Size"


# read & organize files -----
res.alt <- map_dfr(
  readRDS(file_power_all),
  ~as_tibble(.x, rownames = "model") %>%
    pivot_longer(!model, names_to = NULL, values_to = "power") %>%
    mutate("model" = str_extract(model, "\\d+.*\\d+$")),
  .id = "method"
)
res.alt$method <- factor(res.alt$method, c("PCO", "PC1", "minp"), c("Trans-PCO", "PC1", "MinP"))



# boxplot plot & jittered points -----
ggplot(data = res.alt, aes(x = model, y = power, color = method)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  ggbeeswarm::geom_quasirandom(dodge.width = 1, size = 0.5, alpha = 0.5, shape = 16) +
  stat_summary(
    aes(group = method),
    geom = "pointrange",
    fun.data = "mean_cl_normal",
    fun.args = list(conf.int = .95),
    position = position_dodge(width = 1),
    color = "green"
  ) +
  labs(x = xlab_name, y = "Power", color = "Method") +
  scale_colour_manual(
    breaks = c("Trans-PCO", "PC1", "MinP"),
    values = c("Trans-PCO" = "#85192d", "PC1" = "#1d349a", "MinP" = "#e89c31"),
    guide = guide_legend(override.aes = list(size = 0.3))
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
    
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    #legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size= unit(0.5, "cm"),
    
    axis.line = element_line(colour="black"),
    #axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text=element_text(colour = "black", size=12),
    axis.title.y = element_text(angle=90,vjust =2, size=14),
    axis.title.x = element_text(vjust = -0.2, size=14),
    
    plot.margin=unit(c(10,5,5,5),"mm")
  ) +
  coord_cartesian(xlim = c(1.2, n_distinct(res.alt$model)-0.2))

# save figure -----
ggsave(paste0(file_power_all, "3.pdf"), height = 3, width = 4.5)

