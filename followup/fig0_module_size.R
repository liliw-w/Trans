######################################################
#################### Plot #modules v.s. module size ####################
######################################################
rm(list = ls())
library(data.table)
library(tidyverse)

file_coexp_module <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds"
coexp_module <- readRDS(file_coexp_module)$moduleLabels


### modules and their sizes
df_module <- tibble("gene" = names(coexp_module),
                       "module" = coexp_module)
df_module_size <- df_module %>% group_by(module) %>% summarise(module_size = n())

### plot histogram of number of modules v.s. module size
fig_modules_size <- ggplot(df_module_size) +
  geom_histogram(aes(x = module_size, fill = after_stat(count)),
                 binwidth = 10,
                 color = "white") +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(8, "Blues")[8:4],
                       na.value = "red") +
  labs(x = "Module Size", y = "Number of Modules") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        axis.line = element_line(colour="black", size = 0.7),
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.text.x = element_text(hjust=1, vjust = 1, size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(vjust = -0.2, size=14),
        axis.title.y = element_text(angle=90,vjust =2, size=14)
  )

### save figure and plotting object for further editing
saveRDS(fig_modules_size, 'fig0a_module_size.rds')

ggsave('fig0a_module_size.png', fig_modules_size, width = 6, height = 5)

