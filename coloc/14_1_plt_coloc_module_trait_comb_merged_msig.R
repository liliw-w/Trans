############################################################
########## Plot merged coloc regions for each pair of (module, trait) ##########
########## for all types of traits: blood, immune, other ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
dir_coloc_gwas <- list.files(
  "/scratch/midway2/liliw1/coloc_MSigDB",
  pattern = paste0("^ukbb_continuous_\\d+"),
  full.names = TRUE
)
file_list_resColoc <- sapply(
  dir_coloc_gwas,
  list.files,
  pattern = "^coloc_reg_w_merged.txt$", full.names = TRUE, recursive = TRUE
)

file_trait_type <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv'


# read files -----
resColoc_all <- lapply(file_list_resColoc, fread, header = TRUE, drop = "Phenocode")
resColoc_all <- bind_rows(resColoc_all[(lapply(resColoc_all, nrow) %>% unlist()) > 0],
                          .id = "trait_type")
trait_type_blood <- fread(file_trait_type, header = TRUE, sep = ",")



# Figure data prep -----

## assign more specific blood traits groups -----
resColoc_all[, "trait_type"] <- resColoc_all %>%
  left_join(trait_type_blood, by = c("trait" = "Trait Abbreviation")) %>%
  pull(`GWAS Group`)



## summarize the number of un-merged and merged coloc regions for (module, trait) -----
resPlot <- resColoc_all %>%
  group_by(trait_type, trait, Module) %>%
  summarise(num_of_regions = n_distinct(Region),
            num_of_regions_merg = n_distinct(merged_region)) %>%
  ungroup()


## set module order by num of coloc regions and change labels -----
file_coexp_module <- '/project2/xuanyao/llw/MODULES/MSigDB/result/coexp.module.rds'
coexp_module <- readRDS(file_coexp_module)

resPlot <- enframe(coexp_module$moduleName, "module_annot", "Module") %>%
  mutate(Module = paste0("module", Module)) %>%
  right_join(resPlot, by = c("Module"))


Module_order = resPlot %>%
  separate(Module, into = c(NA, "module_num"), sep = "module", remove = FALSE, convert = TRUE) %>%
  group_by(Module, module_annot, module_num) %>%
  summarise(n_coloc_regions = sum(num_of_regions),
            n_coloc_regions_merg = sum(num_of_regions_merg),
            n_coloc_trait = n()) %>%
  ungroup() %>%
  arrange(desc(n_coloc_trait), desc(n_coloc_regions_merg), desc(n_coloc_regions))
resPlot$Module = factor(resPlot$Module,
                        levels = Module_order$Module,
                        labels = paste0(Module_order$module_annot))



## set phenotype order within each trait type by num of coloc regions -----
pheno_order <- resPlot %>%
  group_by(trait_type, trait) %>%
  summarise(n_coloc_regions = sum(num_of_regions),
            n_coloc_regions_merg = sum(num_of_regions_merg),
            n_coloc_m = n()) %>%
  ungroup() %>%
  group_by(trait_type) %>%
  arrange(n_coloc_m, n_coloc_regions_merg, n_coloc_regions, .by_group = TRUE) %>%
  ungroup()
resPlot$trait = factor(resPlot$trait,
                       levels = pheno_order$trait)


## assign colocs to trait type -----
resPlot <- mutate(resPlot,
                  "trait_color" = recode(trait_type,
                                         "Platelets" = "#0028a1",
                                         "Red blood cells" = "#85192d",
                                         "White blood cells" = "#e89c31",
                                         "immune" = "#339900",
                                         "other" = "black"))
axis_text_color <- resPlot %>% slice(match(levels(trait), trait))



# tile plot of merged coloc regions for (module, trait) -----
fig_tile <- ggplot(resPlot, aes(Module, trait)) +
  geom_tile(aes(fill = num_of_regions_merg)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(8, "Blues")[3:8],
                       na.value = "red",
                       limits = c(0, max(c(resPlot$num_of_regions, resPlot$num_of_regions_merg)))) +
  labs(y = NULL, x = "Module", fill = "Number of \n Coloc Regions") +
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
        axis.text.x = element_text(colour = "black", size = 6,
                                   angle = 60, hjust = 1, vjust = 1),
        axis.title.y = element_text(vjust = -0.2, size = 12),
        axis.title.x = element_text(vjust = 2, size = 12))

ggsave(
  "ukbb_all/reg_module_trait_merged.pdf",
  fig_tile,
  width = 6, height = 6
)


# tile plot of un-merged coloc regions for (module, trait) -----
fig_tile <- ggplot(resPlot, aes(Module, trait)) +
  geom_tile(aes(fill = num_of_regions)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(8, "Blues")[3:8],
                       na.value = "red",
                       limits = c(0, max(c(resPlot$num_of_regions, resPlot$num_of_regions_merg)))) +
  labs(y = NULL, x = "Module", fill = "Number of \n Coloc Regions") +
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
        axis.text.x = element_text(colour = "black", size = 6,
                                   angle = 60, hjust = 1, vjust = 1),
        axis.title.y = element_text(vjust = -0.2, size = 12),
        axis.title.x = element_text(vjust = 2, size = 12))

ggsave(
  "ukbb_all/reg_module_trait.pdf",
  fig_tile,
  width = 6, height = 6
)
