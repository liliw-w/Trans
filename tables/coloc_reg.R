##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)

# I/O & paras -----
file_list_resColoc <- c(
  list.files(
    "/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits",
    "^.*coloc_reg_w_merged.txt$",
    full.names = TRUE, recursive = TRUE
  ) %>% set_names("blood"),
  list.files(
    "/project2/xuanyao/llw/coloc/immune_traits",
    "^.*coloc_reg_w_merged.txt$",
    full.names = TRUE, recursive = TRUE
  ) %>% set_names("immune"),
  list.files(
    "/project2/xuanyao/llw/coloc/ukbb_coloc_more_traits",
    "^.*coloc_reg_w_merged.txt$",
    full.names = TRUE, recursive = TRUE
  ) %>% set_names("other")
)
file_trait_type <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv'

## output -----
file_coloc_reg_df <- '/project2/xuanyao/llw/coloc/coloc_reg_all.txt'


# read files -----
resColoc_all <- lapply(file_list_resColoc, fread, header = TRUE, drop = "Phenocode")
resColoc_all <- bind_rows(resColoc_all[(lapply(resColoc_all, nrow) %>% unlist()) > 0],
                          .id = "trait_type")
trait_type_blood <- fread(file_trait_type, header = TRUE, sep = ",")



# assign more specific blood traits groups -----
resColoc_all[resColoc_all$trait_type == "blood", "trait_type"] <- resColoc_all %>%
  filter(trait_type == "blood") %>%
  left_join(trait_type_blood, by = c("trait" = "Trait Abbreviation")) %>%
  pull(`GWAS Group`)


# change col names and order -----
resColoc_all %>%
  select(
    trait_type, trait, Module, 
    Region, Pval, merged_region, 
    nsnps, PP.H4.abf, PP.H3.abf, PP.H2.abf, PP.H1.abf, PP.H0.abf
  ) %>%
  separate(Module, c(NA, "module"), sep = "module", convert = TRUE) %>%
  rename(
    gene_module = module,
    trans_region = Region,
    merged_trans_region = merged_region,
    trans_P = Pval
  ) %>%
  fwrite(
    file_coloc_reg_df,
    quote = FALSE, sep = "\t"
  )

