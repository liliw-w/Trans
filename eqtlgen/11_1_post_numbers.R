###############################################################
########### Numbers of eQTLGen signals ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


### paras and I/O
ratio <- 100
file_signal <- paste0('postanalysis/signal_rm_infl_ratio_', ratio, '.txt')


### read data
signal <- fread(file_signal)


### numbers
cat("There are",
    signal %>% nrow(),
    "significant (trans-eQTL, module) pairs. \n",
    
    
    "There are",
    signal %>% distinct(module) %>% nrow(),
    "modules that have at least one trans-eQTL. \n",
    
    
    "There are",
    signal %>%
      distinct(SNP, .keep_all = TRUE) %>% nrow(),
    "eQTLGen SNPs that are significant trans-eQTL for at least one module. \n"
)

