###############################################################
###########  Look at the replication between DGN signals and eQTLGen trans signals ###########
###############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


########## files and parameters, read files ##########
file_eqtlgen_all_snp <- "/project2/xuanyao/llw/eQTLGen/meta.snp.txt.gz"
file_dgn_all_snp <- "/project2/xuanyao/llw/eQTLGen_DGN/DGN.all_snp.txt"

file_eqtlgen_sig <- '/scratch/midway2/liliw1/eQTGen_est_Sigma/postanalysis/LD.prun.in.chr.module.txt'
file_dgn_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/postanalysis/LD.prun.in.chr.module.perm10.fdr10.txt'

file_eqtlgen_sig_orig <- "/project2/xuanyao/llw/eQTLGen_PCO.lambda.01/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz"


eqtlgen_all_snp <- fread(file_eqtlgen_all_snp, header = TRUE) # 10317 all eQTLGen SNPs
dgn_all_snp <- fread(file_dgn_all_snp, header = FALSE) # 6203169 all DGN SNPs
eqtlgen_sig <- fread(file_eqtlgen_sig, header = FALSE)
dgn_sig <- fread(file_dgn_sig, header = FALSE)
eqtlgen_sig_orig <- fread(file_eqtlgen_sig_orig, header = TRUE)


eqtlgen_all_snp <- eqtlgen_all_snp %>%
  unite("SNP_id", SNPChr, SNPPos, sep = ":") %>%
  count(SNP_id) %>% pull(SNP_id)
dgn_all_snp <- dgn_all_snp %>% pull() %>% unique()
eqtlgen_sig <- eqtlgen_sig %>% pull() %>% unique()
dgn_sig <- dgn_sig %>% pull() %>% unique()



### SNP overlaps
cat("There are",
    length(eqtlgen_all_snp),
    "eQTLGen SNPs in total. \n",
    
    "There are",
    length(dgn_all_snp),
    "DGN SNPs in total. \n",
    
    "There are",
    sum(eqtlgen_all_snp %in% dgn_all_snp),
    "eQTLGen SNPs that are also DGN SNPs. \n"
)


### signal overlaps
cat("There are",
    length(eqtlgen_sig),
    "eQTLGen signals. \n",
    
    "There are",
    length(dgn_sig),
    "DGN signals. \n",
    
    "There are",
    sum(eqtlgen_sig %in% dgn_sig),
    "eQTLGen signals that are also DGN signals. \n"
)



