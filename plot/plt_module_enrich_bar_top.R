###############################################################
########### Bar plot enrichment in module genes ###########
###############################################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggrepel)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
fdr_level <- 0.05
module <- 54
file_enrich <- paste0('/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/M', module, '_gprofiler.csv')

## output -----
file_plt_enrich <- paste0('/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/plot/M', module, '_gprofiler_bar.pdf')


# read files -----
enrich <- fread(file_enrich)
thre_sig <- -log10(fdr_level)

View(enrich)

term_group <- c("GO:BP", "KEGG", "REAC", "GO:MF")
enrich <- arrange(enrich, desc(negative_log10_of_adjusted_p_value)) %>%
  filter(source %in% term_group) %>%
  distinct(term_name, .keep_all = TRUE)
n_term_show <- 10 # nrow(enrich)



# plot enrichment of terms in cis genes -----
enrich$term_name <- fct_reorder(
  enrich$term_name, 
  enrich$negative_log10_of_adjusted_p_value, 
  .desc = FALSE
)


slice(enrich, 1:n_term_show) %>%
  ggplot(aes(y = term_name, x = negative_log10_of_adjusted_p_value)) +
  geom_col(
    aes(fill = negative_log10_of_adjusted_p_value),
    show.legend = FALSE, position = position_dodge()
  ) +
  geom_vline(xintercept = thre_sig, linetype = "dashed", color = "black") +
  labs(
    y = NULL, 
    x = NULL, 
    color = NULL,
    shape = NULL
  ) +
  scale_y_discrete(label = function(x) str_wrap(x, width = 80)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_gradient(low = "#7f99b2", high = "#003366", guide = NULL) +
  theme_my_pub() +
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 7)
  )


# print out key message or write out -----
ggsave(
  file_plt_enrich,
  width = 3, height = 1.5
)

