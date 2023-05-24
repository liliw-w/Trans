###############################################################
########### plot partitioned h2 enrichment ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
module <- 4

file_h2_enrich <- paste0('/project2/xuanyao/llw/ldsc/h2_enrich_comb/M', module, '_blood_traits.results')
file_num_coloc <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data_enrich/all.traits.enrich.module_QTL_sig.txt'

file_fig <- paste0("/project2/xuanyao/llw/ldsc/plots/", basename(file_h2_enrich), ".pdf")


# read data -----
h2_enrich <- fread(file_h2_enrich, header = TRUE, sep = "\t")
num_coloc <- fread(file_num_coloc, header = TRUE)


# add number of coloc regions of (module, trait) -----
num_coloc <- num_coloc %>%
  group_by(Phenocode, trait, Module) %>%
  summarise(num_of_regions = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Module, values_from = num_of_regions, values_fill = 0)
if(! paste0("module", module) %in% colnames(num_coloc)) num_coloc[[paste0("module", module)]] <- 0


h2_enrich <- num_coloc %>%
  select(Phenocode, paste0("module", !!module)) %>%
  rename("trait_id" = Phenocode, "num_coloc" = paste0("module", !!module)) %>%
  right_join(h2_enrich, by = c("trait_id") )


# plot -----
h2_enrich <- arrange(h2_enrich, `GWAS Group`, desc(trait_id))
h2_enrich$trait_id <- factor(h2_enrich$trait_id,
                             levels = h2_enrich$trait_id,
                             labels = paste(h2_enrich$`Trait Abbreviation`))

# 95% CI
xlow <- min(h2_enrich$Enrichment - 1.96*h2_enrich$Enrichment_std_error)
xupp <- max(h2_enrich$Enrichment + 1.96*h2_enrich$Enrichment_std_error) + 1

# title
facet_lab <- paste0("Module ", module)
names(facet_lab) <- unique(h2_enrich$module)

# order traits by enrichment score
h2_enrich$trait_id <- factor(h2_enrich$trait_id,
                             levels = h2_enrich %>%
                               group_by(`GWAS Group`) %>%
                               arrange(desc(Enrichment), .by_group = TRUE) %>%
                               ungroup() %>%
                               pull(trait_id))
# plot error bar
base_fig <- ggplot(h2_enrich,
                   aes(x = Enrichment,
                       y = trait_id,
                       color = `GWAS Group`)) +
  facet_wrap(~module, labeller = labeller(module = facet_lab)) +
  geom_point(size = 2) +
  geom_linerange(aes(xmin = Enrichment - 1.96*`Enrichment_std_error`, xmax = Enrichment + 1.96*`Enrichment_std_error`),
                 size = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#007300") +
  labs(y = NULL, color = "Type", shape = NULL)
if(sum(h2_enrich$num_coloc > 0)){
  base_fig <- base_fig +
    geom_point(data = filter(h2_enrich, num_coloc > 0),
               aes(x = Enrichment + 1.96*`Enrichment_std_error` + 0.5, y = trait_id, shape = "Coloc"),
               size = 1.5, fill = "#328f32", color = "black")
}

base_fig +
  scale_x_continuous(
    limits = c(xlow, xupp),
    breaks = c(0, seq(-100, 100, by = 2))
  ) +
  scale_colour_manual(
    breaks = c("White blood cells", "Red blood cells", "Platelets"),
    values = c("Platelets" = "#0028a1", "Red blood cells" = "#85192d", "White blood cells" = "#e89c31"),
    guide = guide_legend(label.position = "bottom",
                         label.theme = element_text(angle = -90, size = 10))
  ) +
  scale_shape_manual(
    values = 23,
    guide = guide_legend(label.position = "bottom",
                         label.theme = element_text(angle = -90, size = 10),
                         override.aes = list(size = 2))
  ) +
  theme_my_pub(legend.text.size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    panel.grid.major.y = element_line(linetype = "dotted"),
    
    #axis.line = element_blank(),
    
    legend.position = "right",
    legend.background = element_blank(),
    
    strip.text = element_text(face = "bold", size = 12),
    #strip.background = element_rect(colour = "black", fill = NA, size = 1.5),
    strip.background = element_blank(),
    
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )

ggsave(file_fig,
       height = 4.5, width = 3.5)
