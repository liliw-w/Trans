############################################################
########### Plot #Indep null SNPs v.s. module size ###########
########### under a pre-spefified z-score threshold for being null ###########
############################################################
rm(list = ls())
library(tidyverse)
library(ggforce)


### paras and I/O
thre_p_z <- 1e-4
file_null_SNP <- 'null_SNP/num_nullSNP.rds'

# output
file_plot_nullSNP_module_size <- paste0('plot/plot_nullSNP_module_size_pthre', thre_p_z, '.pdf')
file_plot_nullSNP_module_size_line <- paste0('plot/plot_nullSNP_module_size_pthre', thre_p_z, '_ratio.pdf')
file_ratio_module <- 'null_SNP/ratio_module.txt'


### read data
res_nullSNP <- readRDS(file_null_SNP)


### 1. bar plot
# data preparation
fig_dat <- res_nullSNP %>%
  filter(thre_z == thre_p_z) %>%
  mutate(prop_dim = num_nullSNP_indep/module_size) %>%
  select(num_nullSNP_indep, module, module_size, prop_dim) %>%
  pivot_longer(c(num_nullSNP_indep, module_size), names_to = "type", values_to = "num")

# plot
fig <- ggplot(fig_dat, aes(x = module, y = num, fill = type)) +
  geom_bar(stat = "identity", position = position_identity()) +
  facet_zoom(xy = prop_dim < 3, zoom.size = 0.4) +
  labs(x = "Module", y = "Dimension -\n Num of Genes/Null SNPs", fill = " ") +
  scale_fill_brewer(palette = "Paired", direction = -1)

# add theme
fig +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        legend.position = "top")

# save figure
ggsave(file_plot_nullSNP_module_size, height = 5, width = 10)




### 2. line plot of proportion
# group module into groups of different ratio
ratio_thre <- c(5, 10, 50, 100, 150, Inf)
n_panel <- length(ratio_thre)
panel_title <- c(paste0(c(0, ratio_thre[-c(n_panel, n_panel-1)]),
                        "<Ratio<",
                        ratio_thre[-n_panel]),
                 paste0("Ratio>", ratio_thre[n_panel-1]))


# data preparation
fig_dat_line <- fig_dat %>% distinct(module, .keep_all = TRUE) %>% select(module, prop_dim)
fig_dat_line <- cbind(fig_dat_line, 
                      "group_ratio" = apply(fig_dat_line,
                                            1,
                                            function(x) which(as.numeric(x["prop_dim"]) <= ratio_thre)[1])
)
fig_dat_line$module <- factor(fig_dat_line$module, levels = 1:1000)
fig_dat_line$group_ratio <- factor(fig_dat_line$group_ratio,
                                   levels = 1:n_panel,
                                   labels = panel_title)


# plot
fig <- ggplot(fig_dat_line, aes(x = module, y = prop_dim)) +
  geom_segment(aes(xend = module, yend = 0), color = "grey", size = 0.5) +
  geom_point() +
  facet_wrap(~group_ratio, scales = "free") + #, labeller = "label_parsed"
  labs(x = "Module", y = "Ratio -\n Null SNPs/Module size")
fig + 
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.6, linetype = "solid"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed", color = "#e5e5e5"),
        
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 90, size = 6),
        
        axis.ticks = element_blank(),
        
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "#f5f5f5", colour = "black", size = 0.6, linetype = "solid")
  )

# save figure
ggsave(file_plot_nullSNP_module_size_line, height = 5, width = 10)


# save ratios for each module
fwrite(
  fig_dat_line,
  file_ratio_module,
  quote = FALSE, sep = "\t"
)


