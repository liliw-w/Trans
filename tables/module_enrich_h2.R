##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
module <- 4
file_h2_enrich <- paste0('/project2/xuanyao/llw/ldsc/h2_enrich_comb/M', module, '_blood_traits.results')

## output -----
file_h2_enrich_df <- paste0("/project2/xuanyao/llw/ldsc/plots/", basename(file_h2_enrich), ".txt")


# read data -----
h2_enrich <- fread(file_h2_enrich, header = TRUE, sep = "\t")


# change order and rename cols -----
select(
  h2_enrich,
  `GWAS Group`, trait_id, `Trait Abbreviation`, `GWAS Trait`, 
  `Prop._SNPs`, `Prop._h2`, `Prop._h2_std_error`, Enrichment, Enrichment_std_error, Enrichment_p, Coefficient, Coefficient_std_error, `Coefficient_z-score`
) %>%
  arrange(`GWAS Group`, desc(Enrichment)) %>%
  rename(
    trait_type = `GWAS Group`, phenocode = trait_id, trait = `Trait Abbreviation`, trait_info = `GWAS Trait`
  ) %>%
  fwrite(
    file_h2_enrich_df,
    sep = '\t', quote = FALSE
  )

