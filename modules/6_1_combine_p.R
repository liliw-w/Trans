##############################################
########### Combine small p's across all modules and chrs ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)

# I/O & paras -----
file.p <- list.files(
  '/scratch/midway2/liliw1/msig_eqtlgen/p',
  "p.module\\d+.chr\\d+.Sigma_nullz.rds", 
  full.names = TRUE
)
p_include_thre <- 1e-5

## output -----
file_p_all <- 'FDR/p_all_1e5.rds'


# read small p's across all modules and chrs -----
p.obs <- rbindlist(
  lapply(
    file.p,
    function(x){
      tmp_y = readRDS(x);
      print(x);
      ind = tmp_y$p < p_include_thre
      if(!is.null(tmp_y) & sum(ind) > 0){
        tmp_y = tmp_y[ind, ]
        mutate(tmp_y, "module" = strsplit(basename(x), '.', fixed = T)[[1]][2])
      }
    }
  )
)


# print out key message or write out -----
saveRDS(
  p.obs,
  file_p_all
)


# define signal cutoff


# remove modules and chr's

# write out


