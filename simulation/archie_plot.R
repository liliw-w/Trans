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
    archie_cc_value.rds
    archie_selected_gene.rds
    plt_archie_cc_value.pdf
    plt_archie_selected_gene.pdf
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_cc_alt <- args[1]
file_selected_gene_alt <- args[2]

## output -----
file_fig_cc <- args[3]
file_fig_selected_gene <- args[4]



# read files -----
cc_alt <- readRDS(file_cc_alt)
selected_gene_alt <- readRDS(file_selected_gene_alt)


# figures -----
## distribution of cc value across models -----
cc_alt[, 'archie', , , drop = TRUE] %>%
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
apply(selected_gene_alt, c(3, 4), sum) %>%
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

