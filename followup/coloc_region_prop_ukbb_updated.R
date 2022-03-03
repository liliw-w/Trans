rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)

file_pheno_manifest = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"
dir_coloc = file.path( "~/xuanyao_llw/coloc/ukbb_coloc_blood_traits")
dir_all_trait = "~/xuanyao_llw/coloc/ukbb_coloc_blood_traits"


file_module_QTL_signals = '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
range(module_QTL_signals$p)
range(module_QTL_signals$q)


pheno_manifest = fread(file_pheno_manifest)

########## files and parameters ##########
gwasPhenocode_seq = pheno_manifest$phenocode_uniq


########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 'module_QTL_sig'
nsnpsThre = 5
#gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)

#dir.create(file.path(dir_all_trait, "plot"), recursive = TRUE, showWarnings = FALSE)
#dir.create(file.path(dir_all_trait, "data"), recursive = TRUE, showWarnings = FALSE)


setwd(dir_all_trait)
getwd()
file_plot = paste0("plot/all.traits.coloc.region.summary.pvalThre-", pvalThre, ".png")
file_res_coloc_reg_prop = paste0("data/coloc_region_prop_pvalThre-", pvalThre, ".txt")


########## loop over all traits ##########
########## number of regions & coloc regions, proportion of coloc regions ##########
file_resColoc = list.files(dir_coloc, ".*resColoc.txt.gz", full.names = TRUE, recursive = TRUE)
resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

resColoc_all = resColoc_all %>% mutate('if_module_QTL' = Region %in% module_QTL_signals$signal )

res_coloc_reg_prop = resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegion = n(),
            nRegionColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre),
            nRegionPval = sum(if_module_QTL),
            nRegionPvalColoc = sum(if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)) %>%
  ungroup()
res_coloc_reg_prop = res_coloc_reg_prop %>%
  mutate(propColoc = nRegionColoc/nRegion,
         propPvalColoc = nRegionPvalColoc/nRegionPval) %>%
  arrange(desc(propPvalColoc), desc(propColoc)) 


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
  pivot_longer(c(nRegionColoc, nRegionPvalColoc), names_to = "regionTypeColoc", values_to = "nColoc") %>%
  arrange(desc(propPvalColoc), desc(propColoc))
dat_fig_bar_prop$Phenocode = with(dat_fig_bar_prop,
                                  factor(Phenocode,
                                         levels = unique(Phenocode),
                                         labels = paste(trait[!duplicated(Phenocode)],
                                                        Phenocode[!duplicated(Phenocode)],
                                                        sep = "_"))
)

fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = Phenocode)) +
  geom_bar(aes(y = n, fill = regionType), stat = "identity", position=position_dodge()) +
  geom_bar(aes(y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = "Trait", y = "#Regions") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() + theme(text = element_text(size = 4),
                     axis.text.x = element_text(angle = 60, size = 3),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")
fig_bar_prop

# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop = res_coloc_reg_prop %>% select(c(Phenocode, trait, propColoc, propPvalColoc)) %>%
  arrange(desc(propPvalColoc), desc(propColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")
dat_fig_line_prop$Phenocode = with(dat_fig_line_prop,
                                   factor(Phenocode,
                                          levels = unique(Phenocode),
                                          labels = paste(trait[!duplicated(Phenocode)],
                                                         Phenocode[!duplicated(Phenocode)],
                                                         sep = "_"))
)

fig_line_prop <- ggplot(dat_fig_line_prop, aes(x = Phenocode,
                                               y = proportion, group = Type, color = Type)) +
  geom_line() +
  geom_point() +
  labs(x = "Trait", y = "Coloc Proportion", color = "Region type") +
  scale_colour_manual(values = c("blue", "green")) +
  theme_bw() + theme(text = element_text(size = 4),
                     axis.text.x = element_text(angle = 60, size = 3),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")
fig_line_prop

# combine the plots
fig <- plot_grid(fig_bar_prop, fig_line_prop, labels = c('A', "B"), ncol = 1)

# save plot
ggsave(file_plot, fig, width = 10, height = 7)

