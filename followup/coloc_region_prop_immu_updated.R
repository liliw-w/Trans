rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)


########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 'module_QTL_sig'
nsnpsThre = 5
#gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)

dirPlot = "~/xuanyao_llw/coloc/immune_traits/pmid_all"
#dir.create(dirPlot, showWarnings = FALSE)
file_plot = file.path(dirPlot, paste0("all.traits.coloc.region.summary.pvalThre-", pvalThre, ".png"))
file_res_coloc_reg_prop = file.path(dirPlot, paste0("coloc_region_prop_pvalThre-", pvalThre, ".txt"))

gwas_pmid_seq = c(29892013, 31604244, 24390342, 29083406, 30929738, 26502338, 26192919, 26192919, 26192919, 28067908, 28067908, 28067908)
gwas_label_seq = c("AE", "MS", "RA_GWASmeta_European", "Allergy", "ASTHMA", "sle", "IBD", "CD", "UC", "ibd", "cd", "uc")

dir_gwas = file.path("~/xuanyao_llw/coloc/immune_traits/",
                     paste0("pmid", gwas_pmid_seq, "_", gwas_label_seq))
dir_gwas_data = file.path(dir_gwas, "data")


file_module_QTL_signals = '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
range(module_QTL_signals$p)
range(module_QTL_signals$q)

########## loop over all traits ##########
########## number of regions & coloc regions, proportion of coloc regions ##########
file_resColoc = file.path(dir_gwas_data, "resColoc.txt.gz")
resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

resColoc_all = resColoc_all %>% mutate('if_module_QTL' = Region %in% module_QTL_signals$signal )

res_coloc_reg_prop = resColoc_all %>%
  group_by(Phenocode, trait) %>%
  summarise(nRegion = n(),
            nRegionColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre),
            nRegionPval = sum(if_module_QTL),
            nRegionPvalColoc = sum(if_module_QTL & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)) %>%
  #arrange(desc(nRegionPvalColoc), desc(nRegionColoc), desc(nRegionPval), desc(nRegion)) %>%
  ungroup()
res_coloc_reg_prop = res_coloc_reg_prop %>%
  mutate(propColoc = nRegionColoc/nRegion, propPvalColoc = nRegionPvalColoc/nRegionPval) %>%
  arrange(desc(propPvalColoc), desc(propColoc))


########## save results ##########
fwrite(res_coloc_reg_prop, file_res_coloc_reg_prop, quote = FALSE, sep = "\t")



######################################################
########## 2. visualize the proportions ##########
######################################################

# order the traits based on the number of corresponding regions
#res_coloc_reg_prop = rename(res_coloc_reg_prop, c("Phenocode" = "trait", "trait" = "Phenocode"))


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

fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = Phenocode, y = n, fill = regionType)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_bar(aes(x = Phenocode, y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = NULL, y = "#Regions") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, size = 3),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")
fig_bar_prop

# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop = res_coloc_reg_prop %>%
  select(c(Phenocode, trait, propColoc, propPvalColoc)) %>%
  arrange(desc(propPvalColoc), desc(propColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")
dat_fig_line_prop$Phenocode = with(dat_fig_line_prop,
                                   factor(Phenocode,
                                          levels = unique(Phenocode),
                                          labels = paste(trait[!duplicated(Phenocode)],
                                                         Phenocode[!duplicated(Phenocode)],
                                                         sep = "_"))
)

fig_line_prop <- ggplot(dat_fig_line_prop, aes(x = Phenocode, y = proportion, group = Type, color = Type)) +
  geom_line() +
  geom_point() +
  labs(x = NULL, y = "Coloc Proportion", color = "Region type") +
  scale_colour_manual(values = c("blue", "green")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, size = 3),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")
fig_line_prop

# combine the plots
fig <- plot_grid(fig_bar_prop, fig_line_prop, labels = c('A', "B"), ncol = 1)

# save plot
ggsave(file_plot, fig, width = 7, height = 7)
