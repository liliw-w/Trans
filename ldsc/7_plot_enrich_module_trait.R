rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
file_h2_enrich <- list.files(
  '/project2/xuanyao/llw/ldsc/h2_enrich_par/',
  'M\\d+_all_traits.results',
  full.names = TRUE
)

#file_h2_enrich <- paste0('/project2/xuanyao/llw/ldsc/h2_enrich_par/M', module, '_all_traits.results')
file_num_coloc <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data_enrich/all.traits.enrich.module_QTL_sig.txt'

#file_fig <- paste0("/project2/xuanyao/llw/ldsc/plots/", basename(file_h2_enrich), ".pdf")


# read data -----
h2_enrich <- lapply(file_h2_enrich, fread, header = TRUE, sep = "\t") %>% rbindlist(fill = TRUE)

coloc_fig_order <- readRDS('coloc_m_trait_order.rds')

ord_t <- levels(coloc_fig_order$Phenocode)
ord_m <- levels(coloc_fig_order$Module) %>% as.numeric()

h2_enrich$`Trait Abbreviation` <- factor(
  h2_enrich$`Trait Abbreviation`,
  levels = ord_t
)
h2_enrich$module <- factor(
  h2_enrich$module,
  levels = ord_m
)

saveRDS(h2_enrich, 'h2_enrich.rds')


-log10(h2_enrich$Enrichment_p)

fig_tile <- ggplot(h2_enrich, aes(module, `Trait Abbreviation`)) +
  geom_tile(aes(fill = -log10(Enrichment_p))) +
  scale_fill_gradientn(colors = c("white", RColorBrewer::brewer.pal(8, "Blues")[3:8]),
                       na.value = "red") +
  labs(y = NULL, x = "Module", fill = bquote(h2enrich -log[10]) ) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.text.y = element_text(colour = x_axis_text_color[-(1:6), ]$trait_color, size = 8),
        axis.text.x = element_text(colour = "black", size = 8,
                                   angle = 60, hjust = 1, vjust = 1),
        axis.title.y = element_text(vjust = -0.2, size = 12),
        axis.title.x = element_text(vjust = 2, size = 12))
fig_tile

saveRDS(fig_tile, 'h2_enrich_all.rds')
ggsave("h2_enrich_all.pdf", fig_tile, width = 8.5, height = 5)


#h2_enrich <- fread(file_h2_enrich, header = TRUE, sep = "\t")
num_coloc <- fread(file_num_coloc, header = TRUE)
