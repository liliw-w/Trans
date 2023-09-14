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
    /project2/xuanyao/llw/compare_to_archie/archie_results.xlsx
    /project2/xuanyao/llw/compare_to_archie/null_SNP/p_all_archie.txt.gz
    /project2/xuanyao/llw/compare_to_archie/null_SNP/module_archie.rds
    0.05
    /project2/xuanyao/llw/compare_to_archie/null_SNP/plt_archie_pco.pdf
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_signal_archie <- args[1]
file_p_us <- args[2]
file_module_archie <- args[3]
fdr_level <- as.numeric(args[4])


## output -----
file_fig_archie_pco <- args[5]



# read files -----
signal_archie <- readxl::read_excel(file_signal_archie)
p_us <- data.table::fread(file_p_us, header = TRUE)
module_archie <- readRDS(file_module_archie)


# organize data -----
## add adjusted p by bonforonni correction -----
p_us <- left_join(
  p_us,
  count(p_us, Trait, name = 'n_test'),
  by = c('Trait')
) %>%
  mutate(
    "p_adjust" = ifelse(p * n_test > 1, 1, p * n_test)
  )


## combine with archie signals -----
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



# plot p distribution and archie signals -----
# ## nominal p -----
# ggpubr::ggarrange(
#   filter(p_us, p < fdr_level) %>%
#     ggplot() +
#     facet_wrap(~str_glue('{Trait}_{ARCHIE_Component}'), scales = "free", nrow = 1) +
#     geom_point(aes(x = snp, y = -log10(p)), alpha = 0.2) +
#     geom_point(
#       data = var_archie,
#       aes(x = Variants, y = -log10(p)),
#       color = 'maroon', shape = 4
#     ) +
#     geom_hline(yintercept = -log10(fdr_level), linetype = "dashed") +
#     labs(x = "Variant", y = quote(-Log[10](P))) +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank()
#     ),
#   
#   ggplot(var_archie) +
#     facet_wrap(~str_glue('{Trait}_{ARCHIE_Component}'), scales = "free", nrow = 1) +
#     geom_point(aes(x = Variants, y = -log10(p)), color = 'maroon', shape = 4) +
#     geom_hline(yintercept = -log10(fdr_level), linetype = "dashed") +
#     labs(x = "Variant", y = quote(-Log[10](P))) +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank()
#     ),
#   
#   labels = LETTERS[1:2], nrow = 2, ncol = 1
# )


## adjusted p -----
ggpubr::ggarrange(
  filter(p_us, p_adjust < fdr_level) %>%
    ggplot() +
    facet_wrap(~str_glue('{Trait}_{ARCHIE_Component}'), scales = "free", nrow = 1) +
    geom_point(aes(x = snp, y = -log10(p_adjust)), alpha = 0.2) +
    geom_point(
      data = var_archie,
      aes(x = Variants, y = -log10(p_adjust)),
      color = 'maroon', shape = 4
    ) +
    geom_hline(yintercept = -log10(fdr_level), linetype = "dashed") +
    labs(x = "Variant", y = quote(-Log[10](P[adj]))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  
  ggplot(var_archie) +
    facet_wrap(~str_glue('{Trait}_{ARCHIE_Component}'), scales = "free", nrow = 1) +
    geom_point(aes(x = Variants, y = -log10(p_adjust)), color = 'maroon', shape = 4) +
    geom_hline(yintercept = -log10(fdr_level), linetype = "dashed") +
    labs(x = "Variant", y = quote(-Log[10](P[adj]))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  
  ggpubr::ggarrange(
    ggplot() + theme_void(),
    
    count(var_archie, Trait, ARCHIE_Component, if_sig = p_adjust < fdr_level, name = "n_var") %>%
      ggplot(aes(x = str_glue('{Trait}_{ARCHIE_Component}'), y = n_var, fill = if_sig)) +
      geom_col() +
      geom_text(aes(label = n_var), position = "stack") +
      labs(x = "Archie signal", y = "Number of selected variants", fill = "If PCO signal") +
      scale_fill_brewer(palette = "Paired"),
    
    ggplot() + theme_void(),
    
    labels = c("", "C", ""), nrow = 1, ncol = 3, widths = c(1, 2, 1)
  ),
  
  labels = LETTERS[1:2], nrow = 3, ncol = 1
)


# print out key message or write out -----
ggsave(
  filename = file_fig_archie_pco,
  height = 6, width = 8
)

