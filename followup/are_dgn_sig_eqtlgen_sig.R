########## Look at the replication between DGN signals and eQTLGen trans signals ##########
########## To see the relation between DGN signals and complex traits ##########
########## bcs all 10,317 eQTLGen SNPs considered for trans are trait-associated SNPs ##########

rm(list = ls())
library(data.table)
library(tidyverse)

########## files and parameters, read files ##########
file_eqtlgen_all_snp = "/project2/xuanyao/llw/eQTLGen/meta.snp.txt.gz"
file_dgn_all_snp = "/project2/xuanyao/llw/eQTLGen_DGN/DGN.all_snp.txt"
file_eqtlgen_sig = "/project2/xuanyao/llw/eQTLGen_PCO.lambda.01/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz"
file_qtlColocReg = "/project2/xuanyao/llw/coloc/data/qtlColocReg.txt.gz"

eqtlgen_all_snp = fread(file_eqtlgen_all_snp, header = TRUE) # 10317 all eQTLGen SNPs
dgn_all_snp = fread(file_dgn_all_snp, header = FALSE) # 6203169 all DGN SNPs
eqtlgen_sig = fread(file_eqtlgen_sig, header = TRUE)

## use file from coloc for DGN signals, here defined as those whose p<1e-8, as I don't have permutation based signals yet.
qtlColocReg = fread(file_qtlColocReg, header = TRUE)
dgn_sig = qtlColocReg[qtlColocReg$Pval <= 1e-8, ]

eqtlgen_all_snp$snp_id = with(eqtlgen_all_snp, paste(SNPChr, SNPPos, sep = ":"))


########## use unique snps as signals ##########
eqtlgen_uniq_sig = unique(eqtlgen_sig$SNP) # 3,853 unique eQTLGen trans-eQTLs
dgn_uniq_sig = unique(dgn_sig$rsid) # 1,863 unique DGN trans-eQTLs


########## overlapped snps and signals##########
n_overlap_snp = sum(eqtlgen_all_snp$snp_id %in% dgn_all_snp$V1) # 9056
n_overlap_sig = sum(eqtlgen_uniq_sig %in% dgn_uniq_sig) # 27


########## result print out ##########
cat("In total, eQTLGen has", nrow(eqtlgen_all_snp), "SNPs; \n",
    "In total, DGN has", nrow(dgn_all_snp), "SNPs; \n",
    "They have", n_overlap_snp, "overlapped SNPs. \n \n",

    "eQTLGen has", length(eqtlgen_uniq_sig), "signals; \n",
    "DGN has", length(dgn_uniq_sig), "signals; \n",
    "They have", n_overlap_sig, "overlapped signals. \n")

