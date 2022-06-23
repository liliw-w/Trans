############################################################
############### Plot powers for all scenerios on one plot ###############
############### in multiple facets ###############
############################################################
rm(list = ls())

library(data.table)
library(tidyverse)


modify <- function(input){
  x = input[1]
  change= input[2]
  power_emp = readRDS(x)
  DT = data.table("method" = rep(names(power_emp), each=nrow(power_emp[[1]])),
                  "model" = rownames(do.call(rbind, power_emp)),
                  do.call(rbind, power_emp))
  y = melt(DT, id.vars = c("method", "model"), value.name = "power")[, -"variable"]
  
  y$power = as.numeric(y$power)
  y$model = factor(y$model, levels = unique(y$model))
  y$method = factor(y$method, levels = unique(y$method))
  y$change = factor(change)
  return(y)
}


### power results list
file_power_list <- list(
  c('new_Sigma/power.N.lambda0.1.varb1e-3.K101.rds', "Sample Size"),
  c('new_Sigma/power.caus.lambda0.1.varb1e-3.K101.rds', "Causal Proportion"),
  c('new_Sigma/power.var.lambda0.1.varb1e-3.K101.rds', "Genetic Variance")
)


### read data and combine
res.alt <- lapply(file_power_list, modify) |> rbindlist()


### prepare data for plotting
# add mean and se of power
dat2 <- res.alt %>%
  group_by(method, model, change) %>%
  summarise(m = mean(power),
            se = sd(power)/sqrt(n()) ) %>%
  ungroup()

# set method order and labels
dat2$method <- factor(dat2$method,
                      levels = c("PCO", "minp", "PC1"),
                      labels = c("Trans-PCO", "Minp", "PC1-based"))

# set model labels
dat2$model = factor(dat2$model,
                    levels = levels(dat2$model),
                    labels = sapply(levels(dat2$model), function(x) strsplit(x, "=")[[1]][2] )
)

# set scenerios' order (facet order)
dat2$change <- factor(dat2$change,
                      levels = c("Sample Size", "Causal Proportion", "Genetic Variance"))


### plot
ggplot(dat2, aes(x = model, y = m, color = method)) +
  geom_pointrange(aes(ymin = m-se, ymax = m+se),
                  fatten = 2.5, size = 0.7,
                  position = position_dodge(width = 0.3)) +
  labs(x = "Scenarios", y = "Power", color = "Method") +
  facet_wrap(vars(change), scales = "free_x") +
  scale_colour_manual(
    #breaks = c("trans-PCO", "minp", "PC1"),
    #labels = c("Trans-PCO", "Minp", "PC1-based"),
    values = #c("trans-PCO" = "#074e67", "PC1" = "#67074e", "minp" = "#dd9933")
      c("Trans-PCO" = "#85192d", "Minp" = "#e89c31", "PC1-based" = "#1d349a"),
    guide = guide_legend(override.aes = list(size = 0.3))
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.y = element_line(linetype = "dashed", size = 0.8),
    
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", linetype = "dashed"),
    legend.key.size= unit(0.5, "cm"),
    
    axis.line = element_line(colour="black"),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text=element_text(colour = "black", size=14),
    axis.title.y = element_text(angle=90,vjust =2, size=16),
    axis.title.x = element_text(vjust = -0.2, size=16),
    
    plot.margin=unit(c(10,5,5,5),"mm")
  )


ggsave("updated_power_plot_se.pdf", height = 4, width = 9)
