######################################################
#################### Plot module v.s. trait v.s. #coloc regions ####################
#################### Combine ukbb blood traits and immune traits in one figure ####################
######################################################
rm(list = ls())
library(data.table)
library(tidyverse)

pvalThre = 'module_QTL_sig'

file_resEnrich_blood = paste0("/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data_enrich/all.traits.enrich.", pvalThre, ".txt")
file_resEnrich_immu = paste0("/project2/xuanyao/llw/coloc/immune_traits/pmid_all/data_enrich/all.traits.enrich.", pvalThre, ".txt")
file_plot = 'ukbb_immun_module.pdf'
file_plot_num = 'ukbb_immun_module.txt'

### read module v.s. chr v.s. regions file for blood traits & immune traits
# blood traits
resEnrich_blood <- fread(file_resEnrich_blood, header = TRUE)
resEnrich_blood <- mutate(resEnrich_blood,
                          Phenocode_w_id = paste(Phenocode, trait, sep = ";") )
# immune traits
resEnrich_immu <- fread(file_resEnrich_immu, header = TRUE)
resEnrich_immu$trait <- resEnrich_immu$Phenocode
resEnrich_immu <- mutate(resEnrich_immu,
                         Phenocode_w_id = Phenocode )
# combine
resEnrich <- rbind(resEnrich_blood, resEnrich_immu)


### collect #coloc regions (dark greens) for each trait
resPlot = resEnrich %>%
  group_by(Phenocode, trait, Phenocode_w_id, Module) %>%
  summarise(num_of_regions = n()) %>%
  ungroup()

### add colors to x-axis text by trait type
x_axis_text_color <- NULL
{
  # blood traits type
  file_trait_type <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv'
  trait_type <- fread(file_trait_type, header = TRUE, sep = ",")
  
  # add immmne trait type
  trait_type <- rbind(trait_type,
                      tibble("GWAS Group" = "immune",
                             'GWAS ID' = unique(resEnrich_immu$Phenocode),
                             'GWAS Trait' = unique(resEnrich_immu$Phenocode),
                             'Trait Abbreviation' = unique(resEnrich_immu$Phenocode))
  )
  
  ### assign color to each trait group
  df_color <- trait_type %>%
    group_by(`GWAS Group`) %>%
    summarise() %>%
    ungroup() %>%
    mutate("trait_color" = c("#339900", "#0028a1", "#85192d", "#e89c31"))
  
  ### assign color to each trait
  trait_type <- trait_type %>%
    left_join(df_color, "GWAS Group") %>%
    mutate( Phenocode_w_id = paste(`GWAS ID`, `Trait Abbreviation`, sep = ";") ) %>%
    arrange(desc(`GWAS Group`))
  
  ### add trait color column to the main data df for plotting
  resPlot <- resPlot %>%
    left_join(trait_type %>% select('GWAS Group', "trait_color", 'Trait Abbreviation'),
              by = c("trait" = "Trait Abbreviation") )
  
  ### get color for each trait
  x_axis_text_color <- resPlot %>%
    group_by(Phenocode, trait, Phenocode_w_id) %>%
    summarise(trait_color = unique(trait_color)) %>%
    ungroup()
}


### add orders of modules for plotting
# order the y-axis labels of modules by their #coloc traits
Module_order = resPlot %>%
  separate(Module, into = c(NA, "module_num"), sep = "module", remove = FALSE, convert = TRUE) %>%
  group_by(Module, module_num) %>%
  summarise(n_coloc_regions = sum(num_of_regions),
            n_coloc_trait = n()) %>%
  ungroup() %>%
  arrange(desc(n_coloc_trait), desc(n_coloc_regions))
# order the modules and change labels
resPlot$Module = factor(resPlot$Module,
                        levels = Module_order$Module,
                        labels = paste0(Module_order$module_num))

# phenotype order
pheno_order <- resPlot %>%
  group_by(`GWAS Group`, Phenocode, trait) %>%
  summarise(n_coloc_regions = sum(num_of_regions),
            n_coloc_m = n()) %>%
  ungroup() %>%
  group_by(`GWAS Group`) %>%
  arrange(n_coloc_m, n_coloc_regions, .by_group = TRUE) %>%
  ungroup()
resPlot$Phenocode = factor(resPlot$Phenocode,
                           levels = pheno_order$Phenocode,
                           labels = pheno_order$trait)

### order x-axis labels of traits by their trait type, if there are various trait types
if(exists("trait_type")){
  trait_type[`GWAS Group` == "immune", ]$Phenocode_w_id <- trait_type[`GWAS Group` == "immune", `Trait Abbreviation`]
  
  ### obtain the color for x-axis text labels, aligned with order of labels
  x_axis_text_color <- x_axis_text_color[match(levels(resPlot$Phenocode), x_axis_text_color$trait), ]
}



########## tile plot trait v.s. module v.s. coloc region ##########
fig_tile <- ggplot(resPlot, aes(Module, Phenocode)) +
  geom_tile(aes(fill = num_of_regions)) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(8, "Blues")[3:8],
                       na.value = "red") +
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
        axis.text.y = element_text(colour = x_axis_text_color$trait_color, size = 8),
        axis.text.x = element_text(colour = "black", size = 8,
                                   angle = 60, hjust = 1, vjust = 1),
        axis.title.y = element_text(vjust = -0.2, size = 12),
        axis.title.x = element_text(vjust = 2, size = 12))
fig_tile

### save figure object and figure
saveRDS(fig_tile, 'fig4a_ukbb_immun.rds')
ggsave(file_plot, fig_tile, width = 8.5, height = 5)


######################################################
#################### 3. Write the plot data into numerical table ####################
######################################################

### row: trait; column: module; cell: coloc regions
resPlot_num = resEnrich %>%
  group_by(Phenocode, trait, Module) %>%
  summarise("Region" = paste(sort(Region), collapse = ';')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Phenocode, trait),
              names_from = Module,
              values_from = Region)

# sort the modules by the number of corresponding coloc traits
ind_n_coloc_trait = order(colSums(!is.na(resPlot_num)), decreasing = TRUE)
resPlot_num = resPlot_num %>% select( ind_n_coloc_trait )

fwrite(resPlot_num, file_plot_num, quote = FALSE, sep = "\t")

