############## 1. Extract specific traits from ukbb traits manifest file ###############
############## Based on resources: Traits from eQTLGen traits; Traits from Mu et al paper; Powerful (and interesting) ukbb GWASs ##############
############## 2. download the sum stat files extracted ##############
############## 3. Visualize the sample sizes ##############

rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)

############## files and parameters
file_pheno_manifest = "/project2/xuanyao/llw/GWAS/UKB_nealelab/phenotype_manifest.tsv"
if_run_extract = FALSE

### download ukbb sum stat files
if_run_wget_download = FALSE
file_wget_download = "/scratch/midway2/liliw1/coloc/wget_download_cmd.txt"
dir_download = "/project2/xuanyao/llw/GWAS/UKB_nealelab/"

### output files
file_pheno_manifest_subset = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"
file_plot = "/scratch/midway2/liliw1/coloc/phenotype_sub_N.png"


## read files
pheno_manifest = fread(file_pheno_manifest)
colnames(pheno_manifest)


############## Select GWAS traits
if(if_run_extract){
  
  ### 1. check trait of interest
  trait_of_interest = "height"
  grep(trait_of_interest, pheno_manifest$description, ignore.case = TRUE)
  View(pheno_manifest[grep(trait_of_interest, pheno_manifest$description, ignore.case = TRUE), ])
  
  
  ### 2. powerful ukbb gwas (i.e. with large cases size of EUR)
  View(pheno_manifest %>%
         arrange(desc(n_cases_EUR)) %>%
         select(c(1:7, n_cases_EUR, n_controls_EUR)) %>%
         filter(!is.na(n_controls_EUR))
  )
}


############## phenocode of selected traits
ukbb_trait_phenocode = c("30690", "30760",
                         "21001", "23104",
                         "50",
                         "335",
                         "495",
                         "714", "714.1",
                         "555", "555.2",
                         "3761",
                         "695.42",
                         "174",
                         "250.2",
                         "332",
                         "290.11",
                         "401",
                         "411.4",
                         "172",
                         "278",
                         "300",
                         "185")
### multiple categorical, which need to specify coding
ukbb_categ_trait_phenocode = c("6152", "6152", "20126")
ukbb_categ_trait_coding = c("8", "9", "0")


############## Extract traits
pheno_manifest_subset = pheno_manifest %>%
  filter(phenocode %in% ukbb_trait_phenocode |
           (phenocode %in% ukbb_categ_trait_phenocode &
              coding %in% ukbb_categ_trait_coding) ) %>%
  arrange(desc(n_cases_EUR))

# modify the phenocode and description for "categorical" traits, to make them unique and clear
ind_catg = pheno_manifest_subset$trait_type == "categorical"

pheno_manifest_subset = pheno_manifest_subset %>%
  mutate(phenocode_uniq = phenocode,
         trait = str_replace_all(description, "\\s+", "-"))

pheno_manifest_subset[ind_catg, c("phenocode_uniq", "trait")] =
  with(pheno_manifest_subset[ind_catg,],
       data.table(paste(phenocode, coding, sep = "_"),
                  str_replace_all(coding_description, "\\s+", "-") )
       )


## write out files
fwrite(pheno_manifest_subset, file_pheno_manifest_subset, quote = FALSE, sep = "\t")

## download ukbb sum stat
if(if_run_wget_download){
  fwrite(data.table(paste(pheno_manifest_subset$wget, "-P", dir_download)),
         file_wget_download,
         quote = FALSE, sep = "\t", col.names = FALSE)
  cmd = paste("bash ukbb_wget_download.sh", file_wget_download)
  system(cmd)
}


############## Visualize the sample sizes (controls & cases) from extracted traits
## data preparation
plot_pheno_manifest_subset = pheno_manifest_subset %>%
  select(c(trait_type, phenocode_uniq, coding,
           description, coding_description,
           n_cases_EUR, n_controls_EUR))

# pivot the data for bar plot
plot_pheno_manifest_subset = plot_pheno_manifest_subset %>%
  pivot_longer(cols = c("n_cases_EUR", "n_controls_EUR"),
               names_to = "n_type",
               values_to = "n") %>%
  mutate(n = replace_na(n, 0))

## plot
RColorBrewer::brewer.pal(10, "Paired") # check colors from a pallete
fig <- ggplot(plot_pheno_manifest_subset, aes(
  x = factor(phenocode_uniq,
             levels = unique(phenocode_uniq),
             labels = description[!duplicated(phenocode_uniq)])
  )) +
  geom_bar(aes(y = n, fill = n_type), stat = "identity", position=position_dodge()) +
  labs(y = "Sample size (EUR)", x = "", fill = "") +
  scale_fill_manual(values = c("#E31A1C", "#A6CEE3")) +
  theme_bw() + theme(text = element_text(size = 7),
                     axis.text.x = element_text(angle = 60, size = 4),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.position = "top")

ggsave(file_plot, fig, width = 8, height = 4)
