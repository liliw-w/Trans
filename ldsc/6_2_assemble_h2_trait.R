###############################################################
########### assemble par h2 file across modules for a trait ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
gwasPhenocode <- '29083406_Allergy'


files_par_h2 <- list.files('/project2/xuanyao/llw/ldsc/h2_enrich_par', paste0(gwasPhenocode, "_M\\d+", "_baseline.results"), full.names = TRUE)
#file_gwas_trait_info <- "/project2/xuanyao/llw/GWAS/72_traits_GWAS.csv"

file_out <- paste0('h2_enrich_par/T_', gwasPhenocode, '_all_modules.results')


# read data -----
par_h2 <- map_dfr(
  files_par_h2,
  ~fread(cmd = paste("sed -n -e 1p -e 99p", .x), sep = "\t")
)
if(!all(par_h2$Category == "L2_1")) stop("Not all extracted rows are from custom annotation.")

#gwas_trait_info <- fread(file_gwas_trait_info, sep = ",", header = TRUE)


# re-arrange data -----
# add annotation name
par_h2 <- par_h2 %>%
  mutate(Category = str_extract(basename(!!files_par_h2), paste0(!!gwasPhenocode, "_M\\d+"))) %>%
  separate(Category, c("trait_id", "trait", "module"), sep = "_", remove = FALSE, convert = TRUE) %>%
  separate(module, c(NA, "module"), sep = "M", convert = TRUE)


# output -----
fwrite(par_h2, file_out, quote = FALSE, sep = "\t")

