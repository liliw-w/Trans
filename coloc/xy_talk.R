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

## output -----
file_plt_prop <- "/project2/xuanyao/llw/coloc/coloc_prop_all.pdf"



# read files -----
res_coloc_reg_prop1 <- fread(file_res_coloc_reg_prop1, header = TRUE)
res_coloc_reg_prop2 <- fread(file_res_coloc_reg_prop2, header = TRUE)
res_coloc_reg_prop3 <- fread(file_res_coloc_reg_prop3, header = TRUE)

res_coloc_reg_prop2$trait <- res_coloc_reg_prop2$Phenocode

res_coloc_reg_prop <- rbind(res_coloc_reg_prop1, res_coloc_reg_prop2, res_coloc_reg_prop3)


df_plt1 <- res_coloc_reg_prop1 %>%
  arrange(propPvalColocMerg) %>%
  pivot_longer(cols = c(nRegionPvalMerg), names_to = "if_merge", values_to = "Num_reg") %>%
  pivot_longer(cols = c(nRegionPvalColocMerg), names_to = "if_merge_coloc", values_to = "Num_reg_coloc") %>%
  pivot_longer(cols = c(propPvalColocMerg), names_to = "if_merge_prop", values_to = "prop") %>%
  mutate(if_merge = replace(if_merge, if_merge == "nRegionPvalMerg", "(trans)"),
         if_merge_coloc = replace(if_merge_coloc, if_merge_coloc == "nRegionPvalColocMerg", "(coloc)"),
         trait_type = paste("Blood", if_merge),
         trait_type_coloc = paste("Blood", if_merge_coloc),
         trait_type_prop = paste("Blood", if_merge_prop))

df_plt2 <- res_coloc_reg_prop2 %>%
  arrange(propPvalColocMerg) %>%
  pivot_longer(cols = c(nRegionPvalMerg), names_to = "if_merge", values_to = "Num_reg") %>%
  pivot_longer(cols = c(nRegionPvalColocMerg), names_to = "if_merge_coloc", values_to = "Num_reg_coloc") %>%
  pivot_longer(cols = c(propPvalColocMerg), names_to = "if_merge_prop", values_to = "prop") %>%
  mutate(if_merge = replace(if_merge, if_merge == "nRegionPvalMerg", "(trans)"),
         if_merge_coloc = replace(if_merge_coloc, if_merge_coloc == "nRegionPvalColocMerg", "(coloc)"),
         trait_type = paste("Autoimmune", if_merge),
         trait_type_coloc = paste("Autoimmune", if_merge_coloc),
         trait_type_prop = paste("Autoimmune", if_merge_prop))

df_plt3 <- res_coloc_reg_prop3 %>%
  arrange(propPvalColocMerg) %>%
  pivot_longer(cols = c(nRegionPvalMerg), names_to = "if_merge", values_to = "Num_reg") %>%
  pivot_longer(cols = c(nRegionPvalColocMerg), names_to = "if_merge_coloc", values_to = "Num_reg_coloc") %>%
  pivot_longer(cols = c(propPvalColocMerg), names_to = "if_merge_prop", values_to = "prop") %>%
  mutate(if_merge = replace(if_merge, if_merge == "nRegionPvalMerg", "(trans)"),
         if_merge_coloc = replace(if_merge_coloc, if_merge_coloc == "nRegionPvalColocMerg", "(coloc)"),
         trait_type = paste("Other", if_merge),
         trait_type_coloc = paste("Other", if_merge_coloc),
         trait_type_prop = paste("Other", if_merge_prop))

df_plt <- rbind(df_plt3, df_plt2, df_plt1)


# add trait order
df_plt$trait <- fct_inorder(factor(df_plt$trait))


x_lim <- max(unique(df_plt$Num_reg))
x_nudge <- 2


base_plt <- ggplot(df_plt, aes(y = trait)) +
  geom_bar(aes(x = Num_reg, fill = trait_type), stat = "identity", position = position_dodge(width = 1)) +
  geom_bar(aes(x = Num_reg_coloc, fill = trait_type_coloc), stat = "identity", position = position_dodge(width = 1)) +
  #geom_point(aes(x = prop*y_lim, color = trait_type_prop),
  #           #position = position_dodge(width = 1),
  #           show.legend = FALSE,
  #           shape = 18, size = 2,
  #           position = position_nudge(x = x_nudge)) +
  #geom_line(aes(x = prop*x_lim, group = trait_type_prop, color = trait_type_prop),
  #          #position = position_dodge(width = 1),
  #          show.legend = FALSE,
  #          position = position_nudge(x = x_nudge)) +
  labs(y = NULL, x = "Number of regions", fill = "Region", color = NULL)

base_plt +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  scale_x_continuous(limits = c(0, x_lim),
                     sec.axis = sec_axis(~./x_lim, name = "Coloc Proportion")) +
  #scale_fill_brewer(palette = "Paired", direction = -1)
  scale_fill_manual(
    breaks = c(
      "Blood (trans)", "Blood (coloc)",
      "Autoimmune (trans)", "Autoimmune (coloc)", 
      "Other (trans)", "Other (coloc)"
    ),
    values = c("Blood (trans)" = "#f7eff7", "Autoimmune (trans)" = "#dbe8d5", "Other (trans)" = "#e5e5e5",
               "Blood (coloc)" = "#730073", "Autoimmune (coloc)" = "#008000", "Other (coloc)" = "black")
  ) +
  scale_color_manual(
    breaks = c("Blood propPvalColocMerg", "Autoimmune propPvalColocMerg", "Other propPvalColocMerg"),
    values = c("Blood propPvalColocMerg" = "blue",
               "Autoimmune propPvalColocMerg" = "#18e518",
               "Other propPvalColocMerg" = "#ff0000")
  ) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black")
  )

ggsave(
  file_plt_prop,
  width = 6, height = 6
)

