###############################################################
########### plt histogram and qqplot of all p values ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source("~/Trans/followup/theme_my_pub.R")
source('~/Trans/followup/plt_qq_vector.R')


file_p_all <- 'p/p.module_all.Sigma_nullz.rds'

### read data
p_all <- readRDS(file_p_all)


### 1. plot histogram of all p values
fig <- ggplot(p_all, aes(p)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = "P-value", y = "Number of SNPs" )
fig +
  theme_my_pub() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"))

ggsave("plot/p_hist.pdf", width = 6, height = 4)


### 2. qqplot of all p values
fig <- qqplot(p_all$p)
fig

ggsave("plot/p_qqplot.pdf", width = 4, height = 4, useDingbats = TRUE)
