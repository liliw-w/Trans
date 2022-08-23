############################################################
########## Plot h2 enrich for each pair of (module, trait) ##########
########## for all types of traits: blood, immune, other ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
file_h2_enrich <- c(
  list.files(
    '/project2/xuanyao/llw/ldsc/h2_enrich_comb',
    '^M\\d+_blood_traits.results$',
    full.names = TRUE
  ),
  
  list.files(
    '/project2/xuanyao/llw/ldsc/h2_enrich_comb',
    '^T.*_all_modules.results$',
    full.names = TRUE
  )
)

## to keep trait and module orders same as the orders in coloc figure
file_coloc_fig_order <- '/project2/xuanyao/llw/ldsc/plots/coloc_m_trait_order.rds'

Nmodule <- 166


# read files -----
h2_enrich <- lapply(file_h2_enrich,
                    fread,
                    header = TRUE, sep = "\t",
                    select = c("GWAS Group", "trait_id", "Trait Abbreviation", "module", "Enrichment_p")) %>%
  rbindlist(fill = TRUE)

coloc_fig_order <- readRDS(file_coloc_fig_order)


# data wrangling -----
## to make coloc trait order have the same names as ldsc traits -----
coloc_fig_order$trait <- fct_recode(coloc_fig_order$trait,
           BMI = "Body-mass-index-(BMI)")

## assign coloc orders to ldsc module & trait orders -----
ord_t <- levels(coloc_fig_order$trait)
ord_m <- levels(coloc_fig_order$Module) %>% as.numeric()

h2_enrich$`Trait Abbreviation` <- factor(
  h2_enrich$`Trait Abbreviation`,
  levels = ord_t[ord_t %in% h2_enrich$`Trait Abbreviation`]
)
h2_enrich$module <- factor(
  h2_enrich$module,
  levels = c(ord_m, setdiff(1:Nmodule, ord_m))
)


## assign colors for different trait types -----
axis_text_color <- coloc_fig_order %>%
  slice(match(levels(h2_enrich$`Trait Abbreviation`), trait))


# plot h2 enrich for (module, trait) pairs of same coloc orders -----
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
        axis.text.y = element_text(colour = axis_text_color$trait_color, size = 8),
        axis.text.x = element_text(colour = "black", size = 8,
                                   angle = 60, hjust = 1, vjust = 1),
        axis.title.y = element_text(vjust = -0.2, size = 12),
        axis.title.x = element_text(vjust = 2, size = 12))

## save figure -----
ggsave("h2_enrich_all.pdf", fig_tile, width = 13, height = 5)
