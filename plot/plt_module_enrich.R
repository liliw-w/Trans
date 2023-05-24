###############################################################
########### Plot enrichment in cis genes ###########
###############################################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggrepel)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
go_annot <- c('GO:0007249', 'GO:0043122', 'GO:0002376', 'KEGG:04064')
fdr_level <- 0.05
file_enrich <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/M66_gprofiler.csv'

## output -----
file_plt_cis_enrich <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/plot/M66_gprofiler.pdf"


# read files -----
enrich <- fread(file_enrich)
thre_sig <- -log10(fdr_level)


# plot enrichment of terms in cis genes -----
ggplot(enrich, aes(x = term_name, y = negative_log10_of_adjusted_p_value)) +
  geom_point(
    aes(color = negative_log10_of_adjusted_p_value, shape = source), 
    size = 4, alpha = 0.9, 
    show.legend = FALSE
  ) +
  geom_hline(yintercept = thre_sig, linetype = "dashed", color = "grey") +
  geom_text_repel(data = subset(enrich, term_id %in% go_annot),
                  aes(label = str_wrap(term_name, width = 45)),
                  size = 3,
                  segment.colour = "black",
                  min.segment.length = 0,
                  max.overlaps = 10,
                  nudge_y = 0.2,
                  box.padding = 1,
                  segment.curvature = -0.1,
                  segment.angle = 30,
                  direction = "y",
                  hjust = "right",
                  segment.linetype = 6,
                  arrow = arrow(length = unit(0.015, "npc"))) +
  labs(
    x = "Enrichment Term", 
    y = quote(-Log[10](~italic(P)[adjust])), 
    color = NULL,
    shape = NULL
  ) +
  scale_colour_gradient(low = "#d2a5a5", high = "#730000", guide = NULL, breaks = seq(1, 4, by = 0.01)) +
  theme_my_pub() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# print out key message or write out -----
ggsave(
  file_plt_cis_enrich,
  width = 5, height = 3
)

