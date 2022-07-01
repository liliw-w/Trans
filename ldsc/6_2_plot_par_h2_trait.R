###############################################################
########### plot partitioned h2 enrichment ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
gwasPhenocode <- '28067908_ibd'
facet_lab <- "Inflammaroty bowel disease"
# Allergic diseases

file_h2_enrich <- paste0('h2_enrich_par/T_', gwasPhenocode, '_all_modules.results')
file_num_coloc <- '/project2/xuanyao/llw/coloc/immune_traits/pmid_all/data_enrich/all.traits.enrich.module_QTL_sig.txt'

file_fig <- paste0("plots/", basename(file_h2_enrich), ".pdf")


# read data -----
h2_enrich <- fread(file_h2_enrich, header = TRUE, sep = "\t")
num_coloc <- fread(file_num_coloc, header = TRUE)
trait <- h2_enrich$trait[1]


# add number of coloc regions of (module, trait) -----
num_coloc <- num_coloc %>%
  group_by(Phenocode, Module) %>%
  summarise(num_of_regions = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Phenocode, values_from = num_of_regions, values_fill = 0)


h2_enrich <- num_coloc %>%
  select(Module, !!trait) %>%
  rename("num_coloc" = !!trait) %>%
  separate(Module, c(NA, "module"), sep = "module", convert = TRUE) %>%
  right_join(h2_enrich, by = "module", )


# plot -----
h2_enrich$module <- factor(h2_enrich$module)

xlow <- min(h2_enrich$Enrichment - 1.96*h2_enrich$Enrichment_std_error)
xupp <- max(h2_enrich$Enrichment + 1.96*h2_enrich$Enrichment_std_error) + 1

names(facet_lab) <- unique(h2_enrich$trait)


base_fig <- ggplot(h2_enrich,
                   aes(x = Enrichment,
                       y = module,
                       color = num_coloc > 0 & !is.na(num_coloc) )) +
  facet_wrap(~trait, labeller = labeller(trait = facet_lab)) +
  geom_point(size = 3) +
  geom_linerange(aes(xmin = Enrichment - 1.96*`Enrichment_std_error`, xmax = Enrichment + 1.96*`Enrichment_std_error`),
                 size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#007300") +
  labs(y = NULL, color = "Type", shape = NULL)


base_fig +
  scale_x_continuous(
    limits = c(xlow, xupp),
    breaks = c(0, seq(-100, 100, by = 2))
  ) +
  scale_colour_manual(
    breaks = c("TRUE", "FALSE"),
    labels = c("Coloc", "No coloc"),
    values = c("TRUE" = "#005900", "FALSE" = "#8c8c8c")
  ) +
  theme_my_pub() +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    panel.grid.major.y = element_line(linetype = "dotted"),
    
    #axis.line = element_blank(),
    
    legend.position = "right",
    legend.background = element_blank(),
    
    strip.text = element_text(face = "bold", size = 14),
    #strip.background = element_rect(colour = "black", fill = NA, size = 1.5),
    strip.background = element_blank()
  )

ggsave(file_fig,
       height = 3, width = 5)
