###############################################################
########### Power plot when low caus, high var ###########
###############################################################
library(tidyverse)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    /scratch/midway3/liliw1/paper1_sim/power_lowCaus_highb_changecaus_varb0.2_N500_K101.rds
    /scratch/midway3/liliw1/paper1_sim/plt_power_lowCaus_highb_changecaus_varb0.2_N500_K101.pdf
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_power_all <- args[1]
file_fig_power <- args[2]


# read & organize files -----
res.alt <- map_dfr(
  readRDS(file_power_all),
  ~bind_cols(.x) %>% 
    setNames(names(.x)) %>%
    pivot_longer(everything(), names_to = "model", values_to = "power") %>%
    mutate("model" = as.numeric(str_extract(model, "\\d+.*\\d+$"))),
  .id = "method"
)
res.alt$method <- factor(res.alt$method, c("PCO", "PC1", "minp"), c("Trans-PCO", "PC1", "MinP"))
res.alt$model <- factor(res.alt$model, levels = sort(unique(res.alt$model)))



# boxplot plot & jittered points -----
title_spec <- str_extract(file_power_all, 'power_.*\\.rds') %>% 
  str_remove("\\.rds") %>%
  str_replace_all("_", " ") %>%
  str_wrap(width = 30, whitespace_only = FALSE)

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
  labs(
    x = "Model", 
    y = "Power", 
    color = "Method", 
    title = title_spec
  ) +
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
    
    axis.line = element_line(colour = "black"),
    #axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text = element_text(colour = "black", size = 12),
    axis.title.y = element_text(angle = 90,vjust = 2, size = 14),
    axis.title.x = element_text(vjust = -0.2, size = 14),
    
    plot.margin=unit(c(10,5,5,5),"mm")
  ) +
  coord_cartesian(xlim = c(1.2, n_distinct(res.alt$model)-0.2))


# save figure -----
ggsave(
  file_fig_power, 
  height = 3, width = 4.5
)

