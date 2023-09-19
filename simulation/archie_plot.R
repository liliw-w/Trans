##############################################
############ archie cc plot ############
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
    TRUE
    /project2/xuanyao/llw/compare_to_archie/archie_null_cc_value_z.rds
    /project2/xuanyao/llw/compare_to_archie/archie_null_selected_gene_z.rds
    /project2/xuanyao/llw/compare_to_archie/plt_archie_null_cc_value_z.pdf
    /project2/xuanyao/llw/compare_to_archie/plt_archie_null_selected_gene_z.pdf
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

if_null <- as.logical(args[1])
file_cc <- args[2]
file_selected_gene <- args[3]

## output -----
file_fig_cc <- args[4]
file_fig_selected_gene <- args[5]



# read files -----
cc <- readRDS(file_cc)
selected_gene <- readRDS(file_selected_gene)


# figures -----
if(if_null){
  ## distribution of cc value across models -----
  cc[,, drop = TRUE] %>%
    enframe(name = NULL, value = "cc_value") %>%
    ggplot() +
    geom_boxplot(
      aes(x = factor("Null"), y = cc_value),
      outlier.shape = 4, outlier.size = 1, outlier.alpha = 0.2
    ) +
    labs(x = NULL)
  
  ggsave(
    file_fig_cc,
    height = 3, width = 2
  )
  
  
  ## distribution of number of selected genes across models -----
  apply(selected_gene, 2, sum) %>%
    enframe(name = NULL, value = "num_gene") %>%
    ggplot() +
    geom_boxplot(
      aes(x = factor("Null"), y = num_gene),
      outlier.shape = 4, outlier.size = 1, outlier.alpha = 0.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "maroon") +
    labs(x = NULL, y = "Number of selected genes")
  
  ggsave(
    file_fig_selected_gene,
    height = 3, width = 2
  )
  
  
}else{
  ## distribution of cc value across models -----
  cc[, 'archie', , , drop = TRUE] %>%
    as_tibble(rownames = "model") %>%
    pivot_longer(
      cols = !model,
      names_to = NULL, values_to = "cc_value"
    ) %>%
    ggplot() +
    geom_boxplot(
      aes(x = model, y = cc_value)
    )
  
  ggsave(
    file_fig_cc,
    height = 3, width = 4
  )
  
  
  ## distribution of number of selected genes across models -----
  apply(selected_gene, c(3, 4), sum) %>%
    t(.) %>%
    as_tibble() %>%
    pivot_longer(
      cols = everything(),
      names_to = "model", values_to = "num_gene"
    ) %>%
    ggplot() +
    geom_boxplot(
      aes(x = model, y = num_gene)
    ) +
    geom_hline(yintercept = 30, linetype = "dashed", color = "maroon") +
    labs(x = "Model", y = "Number of selected genes")
  
  ggsave(
    file_fig_selected_gene,
    height = 3, width = 4
  )
  
}

