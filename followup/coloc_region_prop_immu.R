rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)


########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 1e-8
nsnpsThre = 5
#gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)

dirPlot = "/scratch/midway2/liliw1/coloc/pmid_all"
file_plot = "all.traits.coloc.region.summary.png"
file_res_coloc_reg_prop = file.path(dirPlot, "coloc_region_prop.txt")

gwas_pmid_seq = c(29892013, 31604244, 24390342, 29083406, 30929738, 26502338, 26192919, 26192919, 26192919, 28067908, 28067908, 28067908)
gwas_label_seq = c("AE", "MS", "RA_GWASmeta_European", "Allergy", "ASTHMA", "sle", "IBD", "CD", "UC", "ibd", "cd", "uc")

dir_gwas = file.path("/scratch/midway2/liliw1/coloc",
                     paste0("pmid", gwas_pmid_seq, "_", gwas_label_seq))
dir_gwas_data = file.path(dir_gwas, "data")



########## loop over all traits ##########
########## number of regions & coloc regions, proportion of coloc regions ##########
file_resColoc = file.path(dir_gwas_data, "resColoc.txt.gz")
resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

res_coloc_reg_prop = resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegion = n(),
            nRegionColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre),
            nRegionPval = sum(Pval <= pvalThre),
            nRegionPvalColoc = sum(Pval <= pvalThre & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)) %>%
  ungroup()
res_coloc_reg_prop = res_coloc_reg_prop %>% mutate(propColoc = nRegionColoc/nRegion, propPvalColoc = nRegionPvalColoc/nRegionPval)


########## save results ##########
fwrite(res_coloc_reg_prop, file_res_coloc_reg_prop, quote = FALSE, sep = "\t")



######################################################
########## 2. visualize the proportions ##########
######################################################

# order the traits based on the number of corresponding regions
res_coloc_reg_prop = rename(res_coloc_reg_prop, c("Phenocode" = "trait", "trait" = "Phenocode"))
res_coloc_reg_prop$trait = with(res_coloc_reg_prop, reorder(trait, -nRegion))

# figure 1: draw bar plot on number of reigons
dat_fig_bar_prop = res_coloc_reg_prop %>%
  pivot_longer(c(nRegion, nRegionPval), names_to = "regionType", values_to = "n") %>%
  pivot_longer(c(nRegionColoc, nRegionPvalColoc), names_to = "regionTypeColoc", values_to = "nColoc")

fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = trait, y = n, fill = regionType)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_bar(aes(x = trait, y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = NULL, y = "#Regions") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, size = 7),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")

# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop = res_coloc_reg_prop %>% select(c(trait, propColoc, propPvalColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")

fig_line_prop <- ggplot(dat_fig_line_prop, aes(x = trait, y = proportion, group = Type, color = Type)) +
  geom_line() +
  geom_point() +
  labs(x = NULL, y = "Coloc Proportion", color = "Region type") +
  scale_colour_manual(values = c("blue", "green")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, size = 7),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")

# combine the plots
fig <- plot_grid(fig_bar_prop, fig_line_prop, labels = c('A', "B"), ncol = 1)

# save plot
ggsave(file_plot, fig, path = dirPlot, width = 7, height = 7)
