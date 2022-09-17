###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Calculate the data-based sigma for a selected module ###########
########### and use it for generating z-scores ###########
###############################################################
module <- as.numeric(snakemake@params[['module']])


# I/O & paras -----
file_ex_var_regressed <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds'
file_coexp_module <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/coexp.module.rds'

## output -----
file_Sigma <- paste0('inflation/Sigma_DGN_module', module, '.rds')


# read files -----
ex_var_regressed <- readRDS(file_ex_var_regressed)
coexp_module <- readRDS(file_coexp_module)


# look at the modules and select one for simulation -----
coexp_module <- coexp_module$moduleLabels

table(coexp_module)
coexp_module[coexp_module==module]
sort(names(coexp_module[coexp_module==module]))


# calculate the correlation matrix using DGN raw data -----
genes_in_module <- sort(names(coexp_module[coexp_module==module]))
Sigma_sim <- cor(ex_var_regressed[, genes_in_module])

Sigma_sim[1:5, 1:5]


# print out key message or write out -----
saveRDS(Sigma_sim, file_Sigma)


