###############################################################
########### Are there inflations in eQTLGen signals? ###########
########### Calculate the data-based sigma for a selected module ###########
########### and use it for generating z-scores ###########
###############################################################

### I/O
module <- as.numeric(snakemake@params[['module']])

file_ex_var_regressed <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds'
file_coexp_module <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/coexp.module.rds'


### read data
ex_var_regressed <- readRDS(file_ex_var_regressed)
coexp_module <- readRDS(file_coexp_module)

### look at the modules and select a desired one
coexp_module <- coexp_module$moduleLabels

table(coexp_module)
coexp_module[coexp_module==module]
sort(names(coexp_module[coexp_module==module]))


# output file depending on K
K <- sum(coexp_module==module)
#file_Sigma <- paste0('inflation/Sigma_DGN_module', module, '_K', K, '.rds')
file_Sigma <- paste0('inflation/Sigma_DGN_module', module, '.rds')



### calculate the correlation matrix
genes_in_module <- sort(names(coexp_module[coexp_module==module]))
Sigma_sim <- cor(ex_var_regressed[, genes_in_module])

Sigma_sim[1:5, 1:5]


### save the selected Sigma
saveRDS(Sigma_sim, file_Sigma)


