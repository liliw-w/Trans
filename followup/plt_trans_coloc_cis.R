##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
plt_dat <- tibble(
  "Trans-" = 179,
  "Cis-e/s" = 60,
  "Cis-e" = 51,
  "Cis-s" = 41
)

## output -----
file_plt_bar <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/trans_cis/trans_coloc_cis.pdf'

# organize data -----
plt_dat <- pivot_longer(
  plt_dat,
  cols = everything(),
  names_to = "type",
  values_to = "n_loci"
) %>%
  arrange(desc(n_loci))
plt_dat$type <- fct_inorder(plt_dat$type)

# bar plot -----
ggplot(plt_dat, aes(x = type, y = n_loci)) +
  geom_bar(stat = "identity", fill = "#d8d8d8", color = "#333333", width = 0.65) +
  geom_text(aes(label = n_loci), nudge_y = 10) +
  labs(x = NULL, y = "Number of loci") +
  theme_my_pub() +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed"),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank()
  )


# print out key message or write out -----
ggsave(
  file_plt_bar,
  width = 3.5, height = 3
)
