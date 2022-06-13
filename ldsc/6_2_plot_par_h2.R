###############################################################
########### plot partitioned h2 enrichment ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/followup/theme_my_pub.R')


# paras and I/O -----
module <- 66

file_h2_enrich <- paste0('h2_enrich_par/M', module, '_all_traits.results')
file_num_coloc <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data_enrich/all.traits.enrich.module_QTL_sig.txt'

file_fig <- paste0("plots/", basename(file_h2_enrich), ".pdf")


# read data -----
h2_enrich <- fread(file_h2_enrich, header = TRUE, sep = "\t")
num_coloc <- fread(file_num_coloc, header = TRUE)


# add number of coloc regions of (module, trait) -----
num_coloc <- num_coloc %>%
  group_by(Phenocode, trait, Module) %>%
  summarise(num_of_regions = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Module, values_from = num_of_regions, values_fill = 0)


h2_enrich <- num_coloc %>%
  select(Phenocode, paste0("module", !!module)) %>%
  rename("trait_id" = Phenocode, "num_coloc" = paste0("module", !!module)) %>%
  right_join(h2_enrich, by = c("trait_id") )


# plot -----
h2_enrich <- arrange(h2_enrich, `GWAS Group`, desc(trait_id))
h2_enrich$trait_id <- factor(h2_enrich$trait_id,
                             levels = h2_enrich$trait_id,
                             labels = paste(h2_enrich$`Trait Abbreviation`))

xlow <- min(h2_enrich$Enrichment - 1.96*h2_enrich$Enrichment_std_error)
xupp <- max(h2_enrich$Enrichment + 1.96*h2_enrich$Enrichment_std_error) + 1

facet_lab <- paste0("Module ", module)
names(facet_lab) <- unique(h2_enrich$module)


base_fig <- ggplot(h2_enrich,
                   aes(x = Enrichment,
                       y = trait_id,
                       color = `GWAS Group`)) +
  facet_wrap(~module, labeller = labeller(module = facet_lab)) +
  geom_point(size = 3) +
  geom_linerange(aes(xmin = Enrichment - 1.96*`Enrichment_std_error`, xmax = Enrichment + 1.96*`Enrichment_std_error`),
                 size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#007300") +
  geom_point(data = filter(h2_enrich, num_coloc > 0),
             aes(x = Enrichment + 1.96*`Enrichment_std_error` + 0.5, y = trait_id, shape = "Coloc"),
             size = 2, fill = "#328f32", color = "black") +
  labs(y = NULL, color = "Type", shape = NULL)


base_fig +
  scale_x_continuous(
    limits = c(xlow, xupp),
    breaks = c(0, seq(-100, 100, by = 2))
  ) +
  scale_colour_manual(
    breaks = c("White blood cells", "Red blood cells", "Platelets"),
    values = c("Platelets" = "#0028a1", "Red blood cells" = "#85192d", "White blood cells" = "#e89c31"),
    guide = guide_legend(label.position = "bottom",
                         label.theme = element_text(angle = -90))
  ) +
  scale_shape_manual(
    values = 23,
    guide = guide_legend(label.position = "bottom",
                         label.theme = element_text(angle = -90),
                         override.aes = list(size = 3))
  ) +
  theme_my_pub() +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    panel.grid.major.y = element_line(linetype = "dotted"),
    
    #axis.line = element_blank(),
    
    legend.position = "right",
    legend.background = element_blank(),
    
    strip.text = element_text(face = "bold", size = 16),
    #strip.background = element_rect(colour = "black", fill = NA, size = 1.5),
    strip.background = element_blank()
    )

ggsave(file_fig,
       height = 7, width = 5.5)
