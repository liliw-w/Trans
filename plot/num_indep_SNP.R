############## plot Figure3B ##############
############## #SNP signals for each (module, chr) pair ##############
##########################################
rm(list = ls())
library(data.table)
library(tidyverse)

file_signal <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'

signal <- fread(file_signal, col.names = c("module_SNP", "p", "q") )


### write unique SNP signals for each (module, chr) pair to file for follow up plink
signal %>%
  separate(module_SNP, c("module", "chr", "pos"), ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_ID", chr, pos, sep = ":", remove = FALSE) %>%
  select(module, chr, SNP_ID) %>%
  group_by(module, chr) %>%
  group_walk(~ fwrite(.x,
                      file.path("/scratch/midway2/liliw1/figures",
                                paste0(.y$module, ".chr", .y$chr, ".txt")),
                      quote = FALSE, col.names = FALSE) )

### Use plink to obtain independent SNP signals for each (module, chr) pair
df_indep_SNP = NULL
file_uniq_list <- list.files(path = "/scratch/midway2/liliw1/figures/fig3/fig3_data/",
                             pattern = "^module.*chr.*.txt$")
for(file_uniq in file_uniq_list){
  file_indep = paste0("indep.", file_uniq)
  command = paste0("bash /home/liliw1/xuanyao_llw/eQTLGen_est_Sigma/indep.SNP.sh ", file_uniq, " ", file_indep)
  system(command)
  df_indep_SNP = rbind(df_indep_SNP, c(
    file_uniq,
    as.numeric(system(paste0("cat ", file_indep, " | wc -l"), intern = TRUE))
  ))
}
fwrite(df_indep_SNP, "num_indep_SNP.txt", quote = FALSE, col.names = FALSE, sep = '\t')

