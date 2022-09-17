##############################################
########### Extract eQTLGen genes analyzed for trans- associations ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_sumstats <- '/project2/xuanyao/llw/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz'


## output -----
file_eqtlgen_genes <- "/project2/xuanyao/llw/eQTLGen/eqtlgen_genes.txt"


# read files -----
sumstats <- fread(
  file_sumstats,
  header = TRUE,
  select = c('Gene', 'GeneSymbol', 'GeneChr', 'GenePos')
)


# extract eqtlgen genes -----
sumstats <- distinct(sumstats, Gene, .keep_all = TRUE) %>%
  filter(
    GeneChr %in% as.character(1:22)
  ) %>%
  mutate(
    GeneChr = as.numeric(GeneChr)
  ) %>%
  arrange(
    GeneChr, GenePos
  )


# print out key message and write out -----
fwrite(
  sumstats,
  file_eqtlgen_genes,
  quote = FALSE, sep = "\t"
)
