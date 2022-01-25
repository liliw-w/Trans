rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)

file_pheno_manifest = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"
dir_coloc = "/scratch/midway2/liliw1/coloc"

pheno_manifest = fread(file_pheno_manifest)

########## files and parameters ##########
gwasPhenocode_seq = pheno_manifest$phenocode_uniq


########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 1e-8
nsnpsThre = 5
#gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)

dir_all_trait = "/scratch/midway2/liliw1/coloc/all_trait/"
dir.create(file.path(dir_all_trait, "plot"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dir_all_trait, "data"), recursive = TRUE, showWarnings = FALSE)

setwd(dir_all_trait)
file_plot = "plot/all.traits.coloc.region.summary.png"
file_res_coloc_reg_prop = "data/coloc_region_prop.txt"


########## loop over all traits ##########
########## number of regions & coloc regions, proportion of coloc regions ##########
file_resColoc = list.files(dir_coloc, "resColoc.txt.gz", full.names = TRUE, recursive = TRUE)
resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

res_coloc_reg_prop = resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegion = n(),
            nRegionColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre),
            nRegionPval = sum(Pval <= pvalThre),
            nRegionPvalColoc = sum(Pval <= pvalThre & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)) %>%
  ungroup()
res_coloc_reg_prop = res_coloc_reg_prop %>%
  mutate(propColoc = nRegionColoc/nRegion,
         propPvalColoc = nRegionPvalColoc/nRegionPval) %>%
  arrang(desc(nRegion))

########## save results ##########
fwrite(res_coloc_reg_prop, file_res_coloc_reg_prop, quote = FALSE, sep = "\t")



######################################################
########## 2. visualize the proportions ##########
######################################################

res_coloc_reg_prop = fread(file_res_coloc_reg_prop)

# order the traits based on the number of corresponding regions
#res_coloc_reg_prop$Phenocode = with(res_coloc_reg_prop, reorder(Phenocode, -nRegion))

# figure 1: draw bar plot on number of reigons
dat_fig_bar_prop = res_coloc_reg_prop %>%
  pivot_longer(c(nRegion, nRegionPval), names_to = "regionType", values_to = "n") %>%
  pivot_longer(c(nRegionColoc, nRegionPvalColoc), names_to = "regionTypeColoc", values_to = "nColoc")



fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = factor(Phenocode,
                                                        levels = unique(Phenocode),
                                                        labels = trait[!duplicated(Phenocode)]),
                                             y = n, fill = regionType)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_bar(aes(x = trait, y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = "Trait", y = "#Regions") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() + theme(text = element_text(size = 4),
                     axis.text.x = element_text(angle = 60, size = 3),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")

# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop = res_coloc_reg_prop %>% select(c(Phenocode, trait, propColoc, propPvalColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")

fig_line_prop <- ggplot(dat_fig_line_prop, aes(x = factor(Phenocode,
                                                          levels = unique(Phenocode),
                                                          labels = trait[!duplicated(Phenocode)]),
                                               y = proportion, group = Type, color = Type)) +
  geom_line() +
  geom_point() +
  labs(x = "Trait", y = "Coloc Proportion", color = "Region type") +
  scale_colour_manual(values = c("blue", "green")) +
  theme_bw() + theme(text = element_text(size = 4),
                     axis.text.x = element_text(angle = 60, size = 3),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")

# combine the plots
fig <- plot_grid(fig_bar_prop, fig_line_prop, labels = c('A', "B"), ncol = 1)

# save plot
ggsave(file_plot, fig, width = 10, height = 7)

