##############################################
########### to visualize eigenvalues ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')

# I/O & paras -----
file.ex.var.regressed <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds'
file.gene.meta <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/gene.meta.txt'
file.coexp.module <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/coexp.module.rds'
module <- 1

## output -----
file_plt_eigen_value <- paste0('/project2/xuanyao/llw/DGN_no_filter_on_mappability/plot/eigen_value_M', module, '.pdf')


# read files -----
datExpr <- readRDS(file.ex.var.regressed)
gene.meta <- read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module <- readRDS(file.coexp.module)


# find genes in the module and module sigmal -----
gene_in_cluster <- data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
gene_w_pos <- merge(gene_in_cluster, gene.meta)
gene_trans <- gene_w_pos[, "gene"]

Sigma <- cor(datExpr[, gene_trans])


# EVD -----
lambdas <- eigen(Sigma)$values


# Point plot of eigenvalues -----
enframe(lambdas) %>%
  ggplot(aes(x = name, y = -log10(value))) +
  geom_point(shape = 1, size = 0.5, alpha = 0.3) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey") +
  labs(x = "PC", y = bquote(-Log[10]~italic((lambda))), title = "Gene module 1, Size 625") +
  theme_my_pub() +
  theme(title = element_text(size = 14))


# print out key message or write out -----
ggsave(
  file = file_plt_eigen_value,
  width = 5, height = 3
)

