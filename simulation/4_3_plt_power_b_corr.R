rm(list = ls())
library(data.table)
library(tidyverse)

file_power_b_corr_03 <- 'new_Sigma/power.N.lambda0.1.varb1e-3.K101_b_corr_caus0.3.rds'
file_power_b_corr_1 <- 'new_Sigma/power.N.lambda0.1.varb1e-3.K101_b_corr_caus1.rds'
file_power_b_indep_1 <- 'new_Sigma/power.N.lambda0.1.varb1e-3.K101_b_indep_caus1.rds'
file_power_b_indep_03 <- '~/xuanyao_llw/simulation_lambda0.1/new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds'

power_b_corr_03 <- readRDS(file_power_b_corr_03)
power_b_corr_1 <- readRDS(file_power_b_corr_1)
power_b_indep_1 <- readRDS(file_power_b_indep_1)
power_b_indep_03 <- readRDS(file_power_b_indep_03)


dat_plt <- rbind(
  imap_dfr(
    power_b_corr_03,
    ~as_tibble(t(.x)) %>%
      pivot_longer(everything(), names_to = "model", values_to = "power"),
    .id = "method"
  ) %>%
    mutate("case" = "b_corr_caus0.3"),
  
  imap_dfr(
    power_b_corr_1,
    ~as_tibble(t(.x)) %>%
      pivot_longer(everything(), names_to = "model", values_to = "power"),
    .id = "method"
  ) %>%
    mutate("case" = "b_corr_caus1"),
  
  imap_dfr(
    power_b_indep_1,
    ~as_tibble(t(.x)) %>%
      pivot_longer(everything(), names_to = "model", values_to = "power"),
    .id = "method"
  ) %>%
    mutate("case" = "b_indep_caus1"),
  
  imap_dfr(
    power_b_indep_03,
    ~as_tibble(t(.x)) %>%
      pivot_longer(everything(), names_to = "model", values_to = "power"),
    .id = "method"
  ) %>%
    mutate("case" = "b_indep_caus0.3")
) %>%
  filter(method!= "PCO") %>%
  filter((method == "minp" & case == "b_indep_caus0.3" & model == "N=800") | method == "PC1")
table(dat_plt$method)


ggplot(data = dat_plt, aes(x = model, y = log10(power), color = case, shape = method)) +
  #geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA, fill = NA) +
  ggbeeswarm::geom_quasirandom(data = subset(dat_plt, power != 0),
                               dodge.width = 1, size = 0.3, alpha = 0.7) +
  labs(x = "Model", y = "Log10(Power)", color = "Setting", shape = "Method")+
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
    
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    #legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size= unit(0.5, "cm"),
    
    axis.line = element_line(colour="black"),
    
    axis.text=element_text(colour = "black", size=12),
    axis.title.y = element_text(angle=90,vjust =2, size=14),
    axis.title.x = element_text(vjust = -0.2, size=14),
    
    plot.margin=unit(c(10,5,5,5),"mm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)),
         shape = guide_legend(override.aes = list(size = 2))
         )


ggsave(filename = "power_b_corr.pdf", width = 6, height = 3)
