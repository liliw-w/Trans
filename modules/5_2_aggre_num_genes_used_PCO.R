##############################################
########### aggregate each module's size of various gene sets ###########
##############################################
library(data.table)
library(tidyverse)


# I/O & paras -----

## output -----
file_num_gene_snp_used <- 'p/num_gene_snp_used_module_all.Sigma_nullz.txt'


# read files -----
num_gene_snp_used <- rbindlist(
  lapply(
    list.files('p', '^num_gene_snp_used_module_\\d+[.]Sigma_nullz[.]txt$', full.names = TRUE),
    fread, header = TRUE)
)


# print out key message and write out -----
fwrite(
  arrange(num_gene_snp_used, module, SNPChr, SNPPos),
  file_num_gene_snp_used,
  sep = "\t",
  quote = FALSE
)

