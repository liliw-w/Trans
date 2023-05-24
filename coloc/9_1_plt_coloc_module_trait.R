######################################################
#################### Plot module v.s. trait v.s. #coloc regions ####################
######################################################
rm(list = ls())
library(data.table)
library(tidyverse)

pvalThre <- 'module_QTL_sig'

file_resEnrich <- paste0("ukbb_all/data_enrich/all.traits.enrich.", pvalThre, ".txt")
file_plot <- paste0("ukbb_all/traits.module.coloc.region.", pvalThre, ".order_by_n_coloc_trait.pdf")
file_plot_num <- paste0("ukbb_all/traits.module.coloc.region.", pvalThre, ".order_by_n_coloc_trait.txt")


########## trait v.s. module v.s. coloc region ##########
resEnrich <- fread(file_resEnrich, header = TRUE)

### collect #coloc regions (dark greens)
resPlot <- resEnrich %>%
  group_by(Phenocode, trait, Module) %>%
  summarise(num_of_regions = n()) %>%
  ungroup() %>%
  mutate(Phenocode_w_id = paste(Phenocode, trait, sep = ";") )

x_axis_text_color <- NULL
{
  ### add colors to x-axis text by trait type
  file_trait_type <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv'
  trait_type <- fread(file_trait_type, header = TRUE, sep = ",")
  
  ### assign color to each trait group
  df_color <- trait_type %>%
    group_by(`GWAS Group`) %>%
    summarise() %>%
    ungroup() %>%
    mutate("trait_color" = c("#000000", "#b41f2e", "grey55"))
  
  ### assign color to each trait
  trait_type <- trait_type %>%
    left_join(df_color, "GWAS Group") %>%
    mutate( Phenocode_w_id = paste(`GWAS ID`, `Trait Abbreviation`, sep = ";") ) %>%
    arrange(desc(`GWAS Group`))
  
  ### add trait color column to the main data df for plotting
  resPlot <- resPlot %>%
    left_join(trait_type %>% select('GWAS Group', "trait_color", 'Trait Abbreviation'),
              by = c("trait" = "Trait Abbreviation") )
  
  x_axis_text_color <- resPlot %>%
    group_by(Phenocode, Phenocode_w_id) %>%
    summarise(trait_color = unique(trait_color)) %>%
    ungroup()
}


### add orders of modules for plotting
# order the y-axis labels of modules by their #coloc traits
Module_order <- resPlot %>%
  separate(Module, into = c(NA, "module_num"), sep = "module", remove = FALSE, convert = TRUE) %>%
  group_by(Module, module_num) %>%
  summarise(n_coloc_regions = sum(num_of_regions),
            n_coloc_trait = n()) %>%
  ungroup() %>%
  arrange(n_coloc_trait, n_coloc_regions )
# change module labels
resPlot$Module <- factor(resPlot$Module,
                        levels = Module_order$Module,
                        labels = paste0("M", Module_order$module_num))

### order x-axis labels of traits by their trait type, if there are various trait types
resPlot$Phenocode <- factor(resPlot$Phenocode)
if(exists("trait_type")){
  resPlot$Phenocode <- factor(resPlot$Phenocode,
                             levels = trait_type$`GWAS ID`,
                             labels =  trait_type$Phenocode_w_id)
  ### obtain the color for x-axis text labels
  x_axis_text_color <- x_axis_text_color[match(levels(resPlot$Phenocode), x_axis_text_color$Phenocode_w_id), ]
}


#low = "#A6CEE3", mid = "#93c3dd", high = "#144c73"
########## tile plot trait v.s. module v.s. coloc region ##########
fig_tile <- ggplot(resPlot, aes(Phenocode, Module)) +
  geom_tile(aes(fill = num_of_regions)) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(8, "Blues")[3:8],
                       na.value = "red") +
  labs(x = NULL, y = "Module", fill = "Number of \n Coloc Regions") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        axis.line = element_line(colour="black"),
        plot.margin = unit(c(10,5,5,5),"mm"),
        axis.text.x = element_text(colour = x_axis_text_color$trait_color,
                                   angle = 60, hjust=1, vjust = 1, size = 9),
        axis.text.y = element_text(colour = "black", size=8),
        axis.title.x = element_text(vjust = -0.2, size=14),
        axis.title.y = element_text(angle=90,vjust =2, size=14) )


### save figure object and figure
ggsave(file_plot, fig_tile, width = 6, height = 5)


######################################################
#################### 3. Write the plot data into numerical table ####################
######################################################

### row: trait; column: module; cell: coloc regions
resPlot_num <- resEnrich %>%
  group_by(Phenocode, trait, Module) %>%
  summarise("Region" = paste(Region, collapse = ';')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Phenocode, trait),
              names_from = Module,
              values_from = Region)

# sort the modules by the number of corresponding coloc traits
ind_n_coloc_trait <- order(colSums(!is.na(resPlot_num)), decreasing = TRUE)
resPlot_num <- resPlot_num %>% select( ind_n_coloc_trait )

fwrite(resPlot_num, file_plot_num, quote = FALSE, sep = "\t")
