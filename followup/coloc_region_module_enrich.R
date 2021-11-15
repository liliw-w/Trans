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

file_qtlColocReg = "/scratch/midway2/liliw1/coloc/data/qtlColocReg.txt.gz"
file_trait_info = "/scratch/midway2/liliw1/coloc/ukbb_blood_traits.csv"
file_coexp_module = "/scratch/midway2/liliw1/DGN_no_filter_on_mappability/result/coexp.module.rds"
dir_data_enrich = "/scratch/midway2/liliw1/coloc/data_enrich"
dirPlot = "/scratch/midway2/liliw1/coloc/plot/"
file_plot = "traits.module.coloc.region.png"

qtlColocReg = fread(file_qtlColocReg, header = TRUE)
trait_info = fread(file_trait_info, sep = ",", header = TRUE)
coexp_module = readRDS(file_coexp_module)$moduleLabels


########## loop over all traits ##########
resEnrich = NULL
for (gwasPhenocode in gwasPhenocode_seq) {

  ########## files and read files ##########
  file_resColoc = paste0("/scratch/midway2/liliw1/coloc/data/pheno", gwasPhenocode, ".resColoc.txt.gz")
  resColoc = fread(file_resColoc, header = TRUE)

  ########## coloc regions corresponding module ##########
  resColoc = qtlColocReg %>% select(c(Signal, Pval)) %>%
    right_join(y = resColoc, by = c("Signal" = "Region")) %>%
    rename("Region" = "Signal")

  tmpresEnrich = resColoc %>% filter(Pval <= pvalThre & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) %>%
    select(c(Region, Pval, PP.H4.abf)) %>%
    separate(col = Region, into = c("Module"), sep = ":", remove = FALSE, convert = TRUE) %>%
    mutate("Phenocode" = gwasPhenocode)


  ########## write results ##########
  resEnrich = rbind(resEnrich, tmpresEnrich)

  cat("Trait", gwasPhenocode, "is done! \n")
}

########## add trait info ##########
resEnrich = trait_info %>% select(c("GWAS ID", "Trait Abbreviation")) %>%
  right_join(y = resEnrich, by = c("GWAS ID" = "Phenocode")) %>%
  rename("Phenocode" = "GWAS ID", "trait" = "Trait Abbreviation")


########## Add genes for module ##########
coexp_module = data.table("Module" = paste0("module", coexp_module), "Gene" = names(coexp_module))
resEnrich$query = paste0(">", apply(resEnrich, 1, function(x) paste(x, collapse = ";") ))
resEnrich$gene = sapply(resEnrich$Module, function(x) paste(coexp_module[coexp_module$Module == x, Gene], collapse = " ") )


########## write out in format of g:profile multiple query ##########
res = resEnrich %>% select(c(Phenocode, trait, query, gene)) %>%
  pivot_longer(c(query, gene), names_to = NULL, values_to = "enrich")
res %>% group_by(Phenocode, trait) %>%
  group_walk(~ fwrite(.x, file.path(dir_data_enrich, paste0("pheno", .y$Phenocode, ".", .y$trait, ".enrich.txt")), quote = FALSE, col.names = FALSE) )


########## trait v.s. module v.s. coloc region ##########
resPlot = resEnrich %>% group_by(Phenocode, trait, Module) %>% summarise(num_of_regions = n())
resPlot$Module = factor(resPlot$Module, levels = paste0("module", 1:length(unique(coexp_module$Module))-1))
resPlot$Phenocode = factor(resPlot$Phenocode, levels = resPlot$Phenocode, labels = paste(resPlot$Phenocode, resPlot$trait, sep = ";") )


########## tile plot trait v.s. module v.s. coloc region ##########
fig_tile <- ggplot(resPlot, aes(factor(Phenocode), Module)) +
  geom_tile(aes(fill = num_of_regions)) +
  scale_fill_gradient(low="grey", high="red") +
  labs(x = "Trait", y = "Module", fill = "#coloc regions") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top",
                     axis.text.x = element_text(angle = 90, size = 7),
                     axis.text.y = element_text(size = 7))

ggsave(file_plot, fig_tile, path = dirPlot, width = 5, height = 7)

