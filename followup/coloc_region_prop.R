rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)


########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 1e-8
nsnpsThre = 5
gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)
dirPlot = "/scratch/midway2/liliw1/coloc/plot/"
file_plot = "all.traits.coloc.region.summary.png"

file_qtlColocReg = "/scratch/midway2/liliw1/coloc/data/qtlColocReg.txt.gz"
file_trait_info = "/scratch/midway2/liliw1/coloc/ukbb_blood_traits.csv"
qtlColocReg = fread(file_qtlColocReg, header = TRUE)
trait_info = fread(file_trait_info, sep = ",", header = TRUE)


########## loop over all traits ##########
resPlot = NULL
for (gwasPhenocode in gwasPhenocode_seq) {

  ########## files and read files ##########
  file_resColoc = paste0("/scratch/midway2/liliw1/coloc/data/pheno", gwasPhenocode, ".resColoc.txt.gz")
  resColoc = fread(file_resColoc, header = TRUE)

  ########## coloc regions proportion ##########
  resColoc = qtlColocReg %>% select(c(Signal, Pval)) %>% right_join(y = resColoc, by = c("Signal" = "Region"))

  nRegion = nrow(resColoc)
  nRegionColoc = sum(resColoc$PP.H4.abf > pp4Thre & resColoc$nsnps >= nsnpsThre)

  nRegionPval = sum(resColoc$Pval <= pvalThre)
  nRegionPvalColoc = sum(resColoc$Pval <= pvalThre & resColoc$PP.H4.abf > pp4Thre & resColoc$nsnps >= nsnpsThre)


  ########## write results ##########
  resPlot = rbind(resPlot, data.table("Phenocode" = gwasPhenocode,
                                      "nRegion" = nRegion, "nRegionColoc" = nRegionColoc,
                                      "nRegionPval" = nRegionPval, "nRegionPvalColoc" = nRegionPvalColoc))

  cat("Trait", gwasPhenocode, "is done! \n")
}

########## colocalized regions proportion ##########
resPlot = resPlot %>% mutate(propColoc = nRegionColoc/nRegion, propPvalColoc = nRegionPvalColoc/nRegionPval)
# add trait info
resPlot = trait_info %>% select(c("GWAS ID", "Trait Abbreviation")) %>%
  right_join(y = resPlot, by = c("GWAS ID" = "Phenocode")) %>%
  rename("Phenocode" = "GWAS ID", "trait" = "Trait Abbreviation")


########## visualize the proportions ##########
# order the traits based on the number of corresponding regions
resPlot$trait = with(resPlot, reorder(trait, -nRegion))

# figure 1: draw bar plot on number of reigons
dat_fig_bar_prop = resPlot %>%
  pivot_longer(c(nRegion, nRegionPval), names_to = "regionType", values_to = "n") %>%
  pivot_longer(c(nRegionColoc, nRegionPvalColoc), names_to = "regionTypeColoc", values_to = "nColoc")

fig_bar_prop <- ggplot(dat_fig_bar_prop, aes(x = trait, y = n, fill = regionType)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_bar(aes(x = trait, y = nColoc, fill = regionTypeColoc), stat = "identity", position=position_dodge()) +
  labs(x = "Trait", y = "#Regions") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 7),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")

# figure 2: draw line plot on the colocalized region proportion
dat_fig_line_prop = resPlot %>% select(c(trait, propColoc, propPvalColoc)) %>%
  pivot_longer(c(propColoc, propPvalColoc), names_to = "Type", values_to = "proportion")

fig_line_prop <- ggplot(dat_fig_line_prop, aes(x = trait, y = proportion, group = Type, color = Type)) +
  geom_line() +
  geom_point() +
  labs(x = "Trait", y = "Coloc Proportion", color = "Region type") +
  scale_colour_manual(values = c("blue", "green")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 8),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")

# combine the plots
fig <- plot_grid(fig_bar_prop, fig_line_prop, labels = c('A', "B"), ncol = 1)

# save plot
ggsave(file_plot, fig, path = dirPlot, width = 10, height = 7)
