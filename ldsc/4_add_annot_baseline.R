###############################################################
########### add new annotations to baseline model ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
prefix <- 'M66_Chr12'

files_annot_baseline <- list.files('data/1000G_Phase3_baselineLD_v2.2_ldscores', "baselineLD.\\d+.annot.gz", full.names = TRUE)
files_ldscore_baseline <- list.files('data/1000G_Phase3_baselineLD_v2.2_ldscores', "baselineLD.\\d+.l2.ldscore.gz", full.names = TRUE)
files_M_5_baseline <- list.files('data/1000G_Phase3_baselineLD_v2.2_ldscores', "baselineLD.\\d+.l2.M_5_50", full.names = TRUE)

files_annot <- list.files('ldsc_annot', paste0(prefix, ".\\d+.annot.gz"), full.names = TRUE)
files_ldscore <- list.files('ldsc_annot', paste0(prefix, ".\\d+.l2.ldscore.gz"), full.names = TRUE)
files_M_5 <- list.files('ldsc_annot', paste0(prefix, ".\\d+.l2.M_5_50"), full.names = TRUE)

dir_comb <- "ldsc_annot_add_baseline/"


# function to combine baseline file and new annotation file -----
comb_base_new <- function(file_base, file_new, dir_comb, type){
  
  # read data
  x = fread(file_base)
  y = fread(file_new)
  
  # check if SNPs from two files align
  if(nrow(x) != nrow(y)) stop("Different SNPs for baseline file and new annotation!")
  
  
  # name of the new annotation as the appended column name
  prefix_new = basename(file_new) %>% str_split(pattern = "[.]") %>% unlist() %>% .[1]
  
  # output file
  file_out = paste0(
    dir_comb,
    prefix_new,
    ".",
    basename(file_base)
  )
  
  # assign annotation name as column name & combine df's
  if(type == "annot"){
    colnames(y) = prefix_new
    res = cbind(x, y)
  }else if(type == "ldscore"){
    colnames(y)[4] = prefix_new
    res = cbind(x, y[, ..prefix_new])
  }else if(type == "M_5"){
    res = cbind(x, y)
    colnames(res) = NULL
  }
  
  # write out combined df to new file
  fwrite(
    res,
    file = file_out,
    quote = FALSE,
    sep = "\t"
  )
}


# 1. append annot file of the baseline file and the new annotation file -----
# iterate over all chr's
walk2(files_annot_baseline, files_annot, comb_base_new, dir_comb, type = "annot")


# 2. append ldscore file of the baseline file and the new annotation file -----
# iterate over all chr's
walk2(files_ldscore_baseline, files_ldscore, comb_base_new, dir_comb, type = "ldscore")


# 3. append M_5 file of the baseline file and the new annotation file -----
# iterate over all chr's
walk2(files_M_5_baseline, files_M_5, comb_base_new, dir_comb, type = "M_5")
