##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')

theme_set(
  theme_my_pub(legend.position = "right")
)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_signal_archie <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/archie/archie_results.xlsx'
file_p_us <- str_glue('null_SNP/p_all_archie.txt.gz')
file_module_archie <- 'null_SNP/module_archie.rds'


## output -----
file_fig_archie_pco <- "null_SNP/plt_archie_pco.pdf"



# read files -----
signal_archie <- readxl::read_excel(file_signal_archie)
p_us <- data.table::fread(file_p_us, header = TRUE)
module_archie <- readRDS(file_module_archie)


# organize data -----
var_archie <- group_by(signal_archie, Trait, ARCHIE_Component) %>%
  select(Variants) %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  
  # add module meta info
  left_join(
    distinct(module_archie, Trait, ARCHIE_Component, module = str_glue('M{module}')),
    by = c('Trait', 'ARCHIE_Component')
  ) %>%
  
  # add trans-pco p
  left_join(
    p_us,
    by = c('Trait', 'ARCHIE_Component', 'module', 'Variants' = 'snp')
  )



# def significance by bonforonni correction, 0.05/#snps, or 0.05/#snps/#gene_sets



# print out key message or write out -----
ggpubr::ggarrange(
  filter(p_us, p < 0.05) %>%
    ggplot() +
    facet_wrap(~str_glue('{Trait}_{ARCHIE_Component}'), scales = "free", nrow = 1) +
    geom_point(aes(x = snp, y = -log10(p)), alpha = 0.2) +
    geom_point(
      data = var_archie,
      aes(x = Variants, y = -log10(p)),
      color = 'maroon', shape = 4
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(x = "Variant", y = quote(-Log[10](P))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  
  ggplot(var_archie) +
    facet_wrap(~str_glue('{Trait}_{ARCHIE_Component}'), scales = "free", nrow = 1) +
    geom_point(aes(x = Variants, y = -log10(p)), color = 'maroon', shape = 4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(x = "Variant", y = quote(-Log[10](P))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  
  labels = LETTERS[1:2], nrow = 2, ncol = 1
)


# print out key message or write out -----
ggsave(
  filename = file_fig_archie_pco,
  height = 4, width = 8
)

