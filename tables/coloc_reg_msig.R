##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)

# I/O & paras -----
file_resColoc_all <- '/project2/xuanyao/llw/MODULES/MSigDB/coloc_MSigDB/ukbb_all/coloc_region_summary_all_merged_annot.txt'
file_trait_type <- '/project2/xuanyao/llw/coloc/ukbb_coloc_blood_traits/ukbb_blood_traits.csv'

## output -----
file_coloc_reg_df <- '/project2/xuanyao/llw/MODULES/MSigDB/coloc_MSigDB/ukbb_all/coloc_region_summary_all_merged_annot_table.txt'

# read files -----
resColoc_all <- fread(file_resColoc_all, header = TRUE)
trait_type_blood <- fread(file_trait_type, header = TRUE, sep = ",")


# change col names and order -----
resColoc_all %>%
  left_join(trait_type_blood, by = c("trait" = "Trait Abbreviation", 'Phenocode' = 'GWAS ID')) %>%
  arrange(`GWAS Group`, module, Chr, Pos) %>%
  select(
    `GWAS Group`, trait, module_annot, module, 
    Region, Pval, merged_region, 
    nsnps, PP.H4.abf, PP.H3.abf, PP.H2.abf, PP.H1.abf, PP.H0.abf
  ) %>%
  rename(
    trait_type = `GWAS Group`,
    module_msig = module_annot,
    gene_module = module,
    trans_region = Region,
    merged_trans_region = merged_region,
    trans_P = Pval
  ) %>%
  fwrite(
    file_coloc_reg_df,
    quote = FALSE, sep = "\t"
  )

