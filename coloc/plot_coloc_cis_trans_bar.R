###############################################################
########### bar plot coloc with cis-e/s numbers ###########
########### use only numbers appear in manuscript ###########
###############################################################
rm(list = ls())


library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
file_upset_inp <- '/project2/xuanyao/llw/coloc/cis/coloc_upset_set.rds'
file_fig <- '/project2/xuanyao/llw/coloc/cis/coloc_trans_cis_bar.pdf'


# read upset input sets -----
listInput <- readRDS(file_upset_inp)


# plot data -----
plt_df <- tibble(
  "Type" = names(listInput),
  "Loci" = lapply(listInput, length) %>% unlist()
)
plt_df$Type <- fct_reorder(plt_df$Type, plt_df$Loci, max, .desc = TRUE)


# plot -----
ggplot(plt_df, aes(x = Type, y = Loci, label = Loci)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(vjust = 0) +
  labs(x = NULL, y = "Number of Loci") +
  scale_x_discrete(
    breaks = c("trans_all", "coloc", "cis_e_cand", "cis_s_cand"),
    labels = c("Tran-All", "Cis-e/s", "Cis-e", "Cis-s")
  ) +
  theme_my_pub(base_size = 12, axis.title.size = 14, axis.text.size = 12) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    axis.line.y = element_blank(),
    axis.ticks = element_blank()
    )

ggsave(file_fig,
       height = 3, width = 3)

