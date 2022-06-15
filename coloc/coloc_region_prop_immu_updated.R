rm(list = ls())
library(data.table)
library(tidyverse)
library(cowplot)

setwd('/project2/xuanyao/llw/coloc/immune_traits/')

########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 'module_QTL_sig'
nsnpsThre = 5
#gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)

file_resColoc = list.files(".", "resColoc.txt.gz", recursive = TRUE)

dirPlot = "~/xuanyao_llw/coloc/immune_traits/pmid_all"
file_plot = file.path(dirPlot, paste0("all.traits.coloc.region.summary.pvalThre-", pvalThre, ".png"))
file_res_coloc_reg_prop = file.path(dirPlot, paste0("coloc_region_prop_pvalThre-", pvalThre, ".txt"))

gwas_pmid_seq = c(29892013, 31604244, 24390342, 29083406, 30929738, 26502338, 26192919, 26192919, 26192919, 28067908, 28067908, 28067908)
gwas_label_seq = c("AE", "MS", "RA_GWASmeta_European", "Allergy", "ASTHMA", "sle", "IBD", "CD", "UC", "ibd", "cd", "uc")


file_module_QTL_signals = '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
range(module_QTL_signals$p)
range(module_QTL_signals$q)


########## loop over all traits ##########
########## number of regions & coloc regions, proportion of coloc regions ##########
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
  mutate(propColoc = nRegionColoc/nRegion, propPvalColoc = nRegionPvalColoc/nRegionPval) %>%
  arrange(desc(propPvalColoc), desc(propColoc))


########## save results ##########
fwrite(res_coloc_reg_prop, file_res_coloc_reg_prop, quote = FALSE, sep = "\t")



######################################################
########## 2. visualize the proportions ##########
######################################################

# order the traits based on the number of corresponding regions
#res_coloc_reg_prop = rename(res_coloc_reg_prop, c("Phenocode" = "trait", "trait" = "Phenocode"))

n_region_pval_all_trait = 255

res_coloc_reg_prop = fread(file_res_coloc_reg_prop)
res_coloc_reg_prop$nRegionPval_all_trait = n_region_pval_all_trait

# figure 1: draw bar plot on number of reigons
dat_fig_bar_prop = res_coloc_reg_prop %>%
  pivot_longer(c(nRegion, nRegionPval), names_to = "regionType", values_to = "n") %>%
  pivot_longer(c(nRegionColoc, nRegionPvalColoc), names_to = "regionTypeColoc", values_to = "nColoc") %>%
  arrange(desc(propPvalColoc), desc(propColoc))
dat_fig_bar_prop$Phenocode = with(dat_fig_bar_prop,
                                  factor(Phenocode,
                                         levels = unique(Phenocode),
                                         labels = sapply(as.character(unique(Phenocode)),
                                                         function(x) strsplit(x, "_", fixed = TRUE)[[1]][1] ))
)

### Don't draw blue regions, only draw green regions, i.e. regions whose lead-SNPs are trans-eQTLs
dat_fig_bar_prop <- dat_fig_bar_prop %>% filter(regionType == "nRegionPval" & regionTypeColoc == "nRegionPvalColoc")


### Draw bar plot (with dual y-axis)
y_lim <- n_region_pval_all_trait
fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = Phenocode)) +
  geom_bar(aes(y = n_region_pval_all_trait), stat = "identity", fill = "#edf5f9") +
  geom_bar(aes(y = n, fill = regionType), stat = "identity", position=position_dodge()) +
  geom_bar(aes(y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = NULL, y = "Regions") +
  scale_fill_brewer(palette = "Paired") + #c("#A6CEE3", "#1F78B4"), c("#B2DF8A", "#33A02C")
  scale_y_continuous(limits = c(0, y_lim),
                     sec.axis = sec_axis(~./y_lim, name = "Coloc Proportion")
  ) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        axis.line = element_line(colour="black"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.text.x = element_text(angle = 60, hjust=1, vjust = 1, size = 10),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.y = element_text(angle=90,vjust =2, size=14),
        axis.title.x = element_text(vjust = -0.2, size=12),
        axis.title.y.right = element_text(angle = 90) )


# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop = res_coloc_reg_prop %>%
  select(c(Phenocode, trait, propColoc, propPvalColoc)) %>%
  arrange(desc(propPvalColoc), desc(propColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")
dat_fig_line_prop$Phenocode = with(dat_fig_line_prop,
                                  factor(Phenocode,
                                         levels = unique(Phenocode),
                                         labels = sapply(as.character(unique(Phenocode)),
                                                         function(x) strsplit(x, "_", fixed = TRUE)[[1]][1] ))
)


### Don't draw blue regions, only draw green regions, i.e. regions whose lead-SNPs are trans-eQTLs
dat_fig_line_prop <- dat_fig_line_prop %>% filter(Type == "propPvalColoc")


### Add the proportion line on the second y-axis of the above bar plot
fig_combined <- fig_bar_prop +
  geom_line(data = dat_fig_line_prop,
            aes(x = Phenocode,
                y = proportion*y_lim, group = Type), color = "blue") +
  geom_point(data = dat_fig_line_prop,
             aes(x = Phenocode,
                 y = proportion*y_lim, group = Type), color = "blue") +
  theme(axis.text.x = element_text(color = "black"))


### save figure object and figure
saveRDS(fig_combined, 'fig4b_immun.rds')

ggsave(file_plot, fig_combined, width = 6, height = 5)
