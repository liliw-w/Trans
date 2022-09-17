###############################################################
########### combine p-values ################
################ across all modules ############
###############################################################
rm(list = ls())
library(data.table)
library(stringr)


### list all files
file_p_list <- list.files("p", "^p.module\\d+\\.Sigma_nullz.rds$", full.names = TRUE)
file_p_list


### combine all files to one single file
p_all <- rbindlist( lapply(file_p_list, function(x) {
  y = readRDS(x)
  #m = as.numeric(strsplit(strsplit(x, ".", fixed = TRUE)[[1]][2], "module")[[1]][2])
  m = as.numeric(str_extract(x, "\\d+"))
  cbind("module" = m, y)
}) )


### write out the combined file
saveRDS(p_all, 'p/p.module_all.Sigma_nullz.rds')




###############################################################
########### check which modules don't have p-value output ################
################ combine p-values across 22 chr for a module ############
###############################################################

### check what modules don't have p values yet
file_p_all_list <- paste0("p/p.module", 1:166, ".Sigma_nullz.rds")
file_p_list <- list.files("p", "^p.module\\d+\\.Sigma_nullz.rds$", full.names = TRUE)

file_p_all_list[!file_p_all_list %in% file_p_list]


### combine chr's for each module
module <- 156
file_p_chr_list <- list.files("p",
                              paste0("^p\\.module", module, "\\.chr\\d+\\.Sigma_nullz.rds$"),
                              full.names = TRUE)
file_p_chr_list
p_module <- rbindlist(lapply(file_p_chr_list, readRDS))

str(p_module)

### save the combined chr's as a whole module
saveRDS(p_module, paste0("p/p.module", module, ".Sigma_nullz.rds"))

