######################################################
#################### Plot #unique SNPs or #indep SNPs ####################
#################### v.s. chr or module ####################
######################################################
rm(list = ls())
library(data.table)
library(tidyverse)


file_signal <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'

signal <- fread(file_signal, col.names = c("module_SNP", "p", "q") )


### write unique SNP signals for each module or chr to file for follow up plink
signal <- signal %>%
  separate(module_SNP, c("module", "chr", "pos"), ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_ID", chr, pos, sep = ":", remove = FALSE) %>%
  select(module, chr, SNP_ID)

signal %>%
  group_by(chr) %>%
  distinct(SNP_ID) %>%
  group_walk(~ fwrite(.x,
                      file.path("fig0_data/",
                                paste0("chr.", .y$chr, ".txt")),
                      quote = FALSE, col.names = FALSE) )

### Use plink to obtain independent SNP signals for each module or chr
df_indep_SNP = NULL
file_uniq_list <- list.files(path = "fig0_data",
                             pattern = "^module.*.txt$",
                             full.names = TRUE)
for(file_uniq in file_uniq_list){
  file_indep = paste0("fig0_data/indep.", basename(file_uniq) )
  command = paste0("bash /home/liliw1/xuanyao_llw/eQTLGen_est_Sigma/indep.SNP.sh ", file_uniq, " ", file_indep)
  system(command)
  df_indep_SNP = rbind(df_indep_SNP, c(
    file_uniq,
    as.numeric(system(paste0("cat ", file_indep, " | wc -l"), intern = TRUE))
  ))
}
fwrite(df_indep_SNP, "fig0_data/num_indep_SNP_module.txt", quote = FALSE, col.names = FALSE, sep = '\t')





