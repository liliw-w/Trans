##############################################
########### Extract eQTLGen SNPs analyzed for trans- associations ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_sumstats <- '/project2/xuanyao/llw/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz'


## output -----
file_eqtlgen_snps <- "/project2/xuanyao/llw/eQTLGen/eqtlgen_snps.txt"


# read files -----
sumstats <- fread(
  file_sumstats,
  header = TRUE,
  select = c('SNP', 'SNPChr', 'SNPPos', 'AssessedAllele', 'OtherAllele')
)


# extract eqtlgen genes -----
sumstats <- distinct(sumstats, SNP, .keep_all = TRUE) %>%
  filter(
    SNPChr %in% as.character(1:22)
  ) %>%
  mutate(
    SNPChr = as.numeric(SNPChr)
  ) %>%
  arrange(
    SNPChr, SNPPos
  )


# print out key message and write out -----
fwrite(
  sumstats,
  file_eqtlgen_snps,
  quote = FALSE, sep = "\t"
)

