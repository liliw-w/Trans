#############################################
### Select a module from DGN gene modules ###
### And use the correlation of this module as the real Sigma for simulations ###
### In this version, the modules are selected from 166 modules of DGN_no_filter_on_mappability ###
#############################################

rm(list = ls())

### files and read objects
file_ex_var_regressed = '~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds'
file_coexp_module = '~/xuanyao_llw/DGN_no_filter_on_mappability/result/coexp.module.rds'
file_Sigma = './new_Sigma/Sigma-new_DGN_module29_K101.rds'

ex_var_regressed = readRDS(file_ex_var_regressed)
coexp_module = readRDS(file_coexp_module)

### look at the modules and select a proper one
coexp_module = coexp_module$moduleLabels
table(coexp_module)

module = 29
K = 101

coexp_module[coexp_module==module]
sort(names(coexp_module[coexp_module==module]))

### calculate the correlation matrix
genes_in_module = sort(names(coexp_module[coexp_module==module]))
Sigma_sim = cor(ex_var_regressed[, genes_in_module])

Sigma_sim[1:5, 1:5]


### save the select Sigma and save it
saveRDS(Sigma_sim, file_Sigma)
