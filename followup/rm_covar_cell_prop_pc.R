###############################################################
########### remove factors of cell type proportions or pcs out of covariates ###########
###############################################################
rm(list = ls())

library(data.table)
library(tidyverse)

covar <- fread('data/covariates.txt')

colnames(covar)[1] <- "fac_name"
fac_cell_prop_pc <- c('X.heme',
                      'Th', 'Tc', 'Tc_act', 'B', 'PC', 'NK', 'NK_act', 'mono', 'DC', 'DC_act', 'neutro',
                      'pc.expr.1', 'pc.expr.2', 'pc.expr.3', 'pc.expr.4', 'pc.expr.5', 'pc.expr.6', 'pc.expr.7', 'pc.expr.8', 'pc.expr.9', 'pc.expr.10')
covar_upd <- filter(covar, ! fac_name %in% !!fac_cell_prop_pc)
colnames(covar_upd)[1] <- ""


fwrite(covar_upd, 'data/covariates_no_cell_prop_pc.txt', quote = FALSE, sep = "\t")
