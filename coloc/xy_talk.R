############################################################
########## Merged coloc regions & proportions for traits ##########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_res_coloc_reg_prop1 <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/data/coloc_region_prop_merged.txt'
file_res_coloc_reg_prop2 <- '/project2/xuanyao/llw/coloc/immune_traits/pmid_all/coloc_region_prop_merged.txt'
file_res_coloc_reg_prop3 <- '/project2/xuanyao/llw/coloc/ukbb_coloc_more_traits/all_trait/data/coloc_region_prop_merged.txt'

#file_plt_prop <- "coloc_region_prop_merged_blood_immu_other.pdf"



# read files -----
res_coloc_reg_prop1 <- fread(file_res_coloc_reg_prop1, header = TRUE)
res_coloc_reg_prop2 <- fread(file_res_coloc_reg_prop2, header = TRUE)
res_coloc_reg_prop3 <- fread(file_res_coloc_reg_prop3, header = TRUE)

res_coloc_reg_prop2$trait <- res_coloc_reg_prop2$Phenocode

res_coloc_reg_prop <- rbind(res_coloc_reg_prop1, res_coloc_reg_prop2, res_coloc_reg_prop3)


df_plt1 <- res_coloc_reg_prop1 %>%
  arrange(desc(propPvalColocMerg)) %>%
  pivot_longer(cols = c(nRegionPvalMerg), names_to = "if_merge", values_to = "Num_reg") %>%
  pivot_longer(cols = c(nRegionPvalColocMerg), names_to = "if_merge_coloc", values_to = "Num_reg_coloc") %>%
  pivot_longer(cols = c(propPvalColocMerg), names_to = "if_merge_prop", values_to = "prop") %>%
  mutate(if_merge = replace(if_merge, if_merge == "nRegionPvalMerg", "(trans)"),
         if_merge_coloc = replace(if_merge_coloc, if_merge_coloc == "nRegionPvalColocMerg", "(coloc)"),
         trait_type = paste("Blood", if_merge),
         trait_type_coloc = paste("Blood", if_merge_coloc),
         trait_type_prop = paste("Blood", if_merge_prop))

df_plt2 <- res_coloc_reg_prop2 %>%
  arrange(desc(propPvalColocMerg)) %>%
  pivot_longer(cols = c(nRegionPvalMerg), names_to = "if_merge", values_to = "Num_reg") %>%
  pivot_longer(cols = c(nRegionPvalColocMerg), names_to = "if_merge_coloc", values_to = "Num_reg_coloc") %>%
  pivot_longer(cols = c(propPvalColocMerg), names_to = "if_merge_prop", values_to = "prop") %>%
  mutate(if_merge = replace(if_merge, if_merge == "nRegionPvalMerg", "(trans)"),
         if_merge_coloc = replace(if_merge_coloc, if_merge_coloc == "nRegionPvalColocMerg", "(coloc)"),
         trait_type = paste("Autoimmune", if_merge),
         trait_type_coloc = paste("Autoimmune", if_merge_coloc),
         trait_type_prop = paste("Autoimmune", if_merge_prop))

df_plt3 <- res_coloc_reg_prop3 %>%
  arrange(desc(propPvalColocMerg)) %>%
  pivot_longer(cols = c(nRegionPvalMerg), names_to = "if_merge", values_to = "Num_reg") %>%
  pivot_longer(cols = c(nRegionPvalColocMerg), names_to = "if_merge_coloc", values_to = "Num_reg_coloc") %>%
  pivot_longer(cols = c(propPvalColocMerg), names_to = "if_merge_prop", values_to = "prop") %>%
  mutate(if_merge = replace(if_merge, if_merge == "nRegionPvalMerg", "(trans)"),
         if_merge_coloc = replace(if_merge_coloc, if_merge_coloc == "nRegionPvalColocMerg", "(coloc)"),
         trait_type = paste("Other", if_merge),
         trait_type_coloc = paste("Other", if_merge_coloc),
         trait_type_prop = paste("Other", if_merge_prop))

df_plt <- rbind(df_plt1, df_plt2, df_plt3)


# add trait order
df_plt$trait <- fct_inorder(factor(df_plt$trait))


y_lim <- max(unique(df_plt$Num_reg))
y_nudge <- 2


base_plt <- ggplot(df_plt, aes(x = trait)) +
  geom_bar(aes(y = Num_reg, fill = trait_type), stat = "identity", position = position_dodge(width = 1)) +
  geom_bar(aes(y = Num_reg_coloc, fill = trait_type_coloc), stat = "identity", position = position_dodge(width = 1)) +
  geom_point(aes(y = prop*y_lim, color = trait_type_prop),
             #position = position_dodge(width = 1),
             show.legend = FALSE,
             shape = 18, size = 2,
             position = position_nudge(y = y_nudge)) +
  geom_line(aes(y = prop*y_lim, group = trait_type_prop, color = trait_type_prop),
            #position = position_dodge(width = 1),
            show.legend = FALSE,
            position = position_nudge(y = y_nudge)) +
  labs(x = NULL, y = "Number of regions", fill = "Region", color = NULL)

base_plt +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(0, y_lim),
                     sec.axis = sec_axis(~./y_lim, name = "Coloc Proportion")) +
  #scale_fill_brewer(palette = "Paired", direction = -1)
  scale_fill_manual(
    breaks = c(
      "Blood (trans)", "Blood (coloc)",
      "Autoimmune (trans)", "Autoimmune (coloc)", 
      "Other (trans)", "Other (coloc)"
    ),
    values = c("Blood (trans)" = "#d9eaf3", "Autoimmune (trans)" = "#e2f3d4", "Other (trans)" = "#fde0e0",
               "Blood (coloc)" = "#186191", "Autoimmune (coloc)" = "#24731f", "Other (coloc)" = "#9e1213")
  ) +
  scale_color_manual(
    breaks = c("Blood propPvalColocMerg", "Autoimmune propPvalColocMerg", "Other propPvalColocMerg"),
    values = c("Blood propPvalColocMerg" = "blue",
               "Autoimmune propPvalColocMerg" = "#18e518",
               "Other propPvalColocMerg" = "#ff0000")
  ) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        #legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        
        axis.line = element_line(colour="black"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        
        #legend.margin = margin(-0.5,0,0,0, unit="cm"),
        
        axis.text.x = element_text(angle = 60, hjust=1, vjust = 1, size = 8, color = "black"),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.y.right = element_text(angle = 90),
        axis.title.x = element_text(size = 12))

ggsave(
  "coloc_prop_xy.pdf",
  width = 9, height = 3
)
