##############################################
########### coloc with cis-e/s in DGN ###########
##############################################
# load packages -----
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_coloc_cise <- '/project2/xuanyao/llw/coloc/cis/cis_e/data/coloc_reg_w_merged.txt'
file_coloc_ciss <- '/project2/xuanyao/llw/coloc/cis/cis_s/data/coloc_reg_w_merged.txt'

## output -----
file_coloc_cis_df <- '/project2/xuanyao/llw/coloc/cis/coloc_cis.txt'


# read files -----
coloc_cise <- fread(file_coloc_cise, header = TRUE)
coloc_ciss <- fread(file_coloc_ciss, header = TRUE)


# change col names and order -----
coloc_cis <- rbind(
  mutate(coloc_cise, cis_type = "cis_e"),
  mutate(coloc_ciss, cis_type = "cis_s")
)
coloc_cis <- separate(coloc_cis, Module, c(NA, "module"), sep = "module", convert = TRUE) %>%
  arrange(module, Chr, Pos) %>%
  select(cis_type, module, Region, merged_region, gene, Pval, 
         nsnps, PP.H4.abf, PP.H3.abf, PP.H2.abf, PP.H1.abf, PP.H0.abf) %>%
  rename(
    gene_module = module,
    trans_region = Region,
    merged_trans_region = merged_region,
    trans_P = Pval
  )

# print out key message or write out -----
fwrite(
  coloc_cis,
  file_coloc_cis_df,
  quote = FALSE, sep = "\t"
)

