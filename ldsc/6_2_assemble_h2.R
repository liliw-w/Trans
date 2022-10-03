###############################################################
########### assemble par h2 file across traits for a module ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# paras and I/O -----
module <- 66

files_par_h2 <- list.files('h2_enrich_par', paste0("^\\d+_M", module, "_baseline.results"), full.names = TRUE)
file_gwasTraitInfo <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv'

## output -----
file_out <- paste0('h2_enrich_par/M', module, '_all_traits.results')


# read data -----
par_h2 <- map_dfr(
  files_par_h2,
  ~fread(cmd = paste("sed -n -e 1p -e 99p", .x), sep = "\t")
)
if(!all(par_h2$Category == "L2_1")) stop("Not all extracted rows are from custom annotation.")

gwasTraitInfo <- fread(file_gwasTraitInfo, sep = ",", header = TRUE)


# re-arrange data -----
# add annotation name & trait info
par_h2 <- par_h2 %>%
  mutate(Category = str_extract(basename(!!files_par_h2), paste0("^\\d+_M", !!module))) %>%
  separate(Category, c("trait_id", "module"), sep = "_", remove = FALSE, convert = TRUE) %>%
  separate(module, c(NA, "module"), sep = "M", convert = TRUE) %>%
  left_join(gwasTraitInfo, by = c("trait_id" = "GWAS ID")) %>%
  relocate(`GWAS Group`, `GWAS Trait`, `Trait Abbreviation`, .after = trait_id)


# output -----
fwrite(par_h2, file_out, quote = FALSE, sep = "\t")

