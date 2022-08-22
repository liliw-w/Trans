############################################################
########## Merged coloc regions & proportions for traits ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_res_coloc_reg_prop <- 'ukbb_coloc_blood_traits/data/coloc_region_prop_merged.txt'

file_plt_prop <- "ukbb_coloc_blood_traits/plot/coloc_region_prop_merged.pdf"


# read files -----
res_coloc_reg_prop <- fread(file_res_coloc_reg_prop, header = TRUE)

df_plt <- res_coloc_reg_prop %>%
  arrange(desc(propPvalColocMerg)) %>%
  pivot_longer(cols = c(nRegionPval, nRegionPvalMerg), names_to = "if_merge", values_to = "Num_reg") %>%
  pivot_longer(cols = c(nRegionPvalColoc, nRegionPvalColocMerg), names_to = "if_merge_coloc", values_to = "Num_reg_coloc") %>%
  pivot_longer(cols = c(propPvalColoc, propPvalColocMerg), names_to = "if_merge_prop", values_to = "prop")


# add trait order
df_plt$trait <- fct_inorder(factor(df_plt$trait))

y_lim <- max(unique(df_plt$Num_reg))

base_plt <- ggplot(df_plt, aes(x = trait)) +
  geom_bar(aes(y = Num_reg, fill = if_merge), stat = "identity", position = position_dodge(width = 1)) +
  geom_bar(aes(y = Num_reg_coloc, fill = if_merge_coloc), stat = "identity", position = position_dodge(width = 1)) +
  geom_point(aes(y = prop*y_lim, color = if_merge_prop), position = position_dodge(width = 1)) +
  geom_line(aes(y = prop*y_lim, group = if_merge_prop, color = if_merge_prop), position = position_dodge(width = 1)) +
  labs(x = NULL, y = "Number of regions", fill = "Region", color = "Proportion")

base_plt +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(0, y_lim),
                     sec.axis = sec_axis(~./y_lim, name = "Coloc Proportion")) +
  scale_fill_manual(
    breaks = c("nRegionPval", "nRegionPvalColoc", "nRegionPvalMerg", "nRegionPvalColocMerg"),
    values = c("nRegionPval" = "#c0dceb", "nRegionPvalColoc" = "#1b6ca2",
               "nRegionPvalMerg" = "#d0ebb8", "nRegionPvalColocMerg" = "#2d9027")
  ) +
  scale_color_manual(
    breaks = c("propPvalColoc", "propPvalColocMerg"),
    values = c("propPvalColoc" = "blue",
               "propPvalColocMerg" = "#18e518")
  ) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        
        axis.line = element_line(colour="black"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        
        axis.text.x = element_text(angle = 60, hjust=1, vjust = 1, size = 8, color = "black"),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.y.right = element_text(angle = 90),
        axis.title.x = element_text(size = 12))

ggsave(
  file_plt_prop,
  width = 8, height = 3
)
