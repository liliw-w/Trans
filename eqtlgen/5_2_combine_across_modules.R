############################################################
########### Combine the results of number of genes used for SNPs ###########
########### across all modules ###########
############################################################
rm(list = ls())
library(data.table)

### files list of each module
file_list <- list.files(path = "p",
                        pattern = "^num_gene_snp_used_module[0-9]*.Sigma_nullz.txt$",
                        full.names = TRUE)

### combine across modules
num_gene_snp_used_all <- lapply(file_list, fread, header = TRUE, sep = "\t", quote = FALSE) %>%
  rbindlist() %>%
  arrange(module)


### save result
fwrite(num_gene_snp_used_all,
       'p/num_gene_snp_used_module_all.Sigma_nullz.txt',
       sep = "\t", quote = FALSE)
