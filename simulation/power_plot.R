###############################################################
########### Power plot of three methods across simulated scenarios ###########
###############################################################
rm(list = ls())

library(tidyverse)


# files & paras -----
# which case to plot
s <- 4


file_power_all_seq <- c(
  'new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds',
  'new_Sigma/power.caus.lambda0.1.varb1e-3.K101.rds',
  'new_Sigma/power.var.lambda0.1.varb1e-3.K101.rds',
  'new_Sigma/power.pc1.rds'
)
xlab_name_seq <- c(
  "Sample Size",
  "Causal Proportion",
  "Genetic Variance",
  "Sample Size"
)

# read & organize files -----
file_power_all <- file_power_all_seq[s]
xlab_name <- xlab_name_seq[s]

res.alt <- map_dfr(
  readRDS(file_power_all),
  ~as_tibble(.x, rownames = "model") %>%
    pivot_longer(!model, names_to = NULL, values_to = "power") %>%
    mutate("model" = str_extract(model, "\\d+.*\\d+$")),
  .id = "method"
)
res.alt$method <- factor(res.alt$method, c("PCO", "PC1", "minp"), c("Trans-PCO", "PC1", "MinP"))


# point range & line plot ----
ggplot(data = res.alt, aes(x = model, y = power, color = method, group = method)) +
  stat_summary(
    geom = "pointrange",
    fun.data = "mean_cl_normal",
    fun.args = list(conf.int = .95),
    position = position_dodge(width = 0.3),
    fatten = 1, size = 0.5
  ) +
  stat_summary(
    geom = "line",
    fun = "mean",
    position = position_dodge(width = 0.3),
    alpha = 0.5, show.legend = FALSE
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
  coord_cartesian(xlim = c(1.2, n_distinct(res.alt$model) - 0.2))

# save figure -----
ggsave(paste0(file_power_all, "3.pdf"), height = 3, width = 4.5)
