##############################################
###########  ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file_eqtlgen <- '/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz'

## output -----
file_eqtlgen_gene <- '/project2/xuanyao/data/eQTLGen/gene_meta_trans.txt.gz'
file_eqtlgen_snp <- '/project2/xuanyao/data/eQTLGen/snp_meta_trans.txt.gz'



# read files -----
eqtlgen <- data.table::fread(
  file_eqtlgen,
  select = c('SNP', 'SNPChr', 'SNPPos', 'Gene', 'GeneSymbol', 'GeneChr',	'GenePos')
)


# organize data -----
distinct(eqtlgen, Gene, GeneSymbol, GeneChr, GenePos) %>%
  data.table::fwrite(
    .,
    file = file_eqtlgen_gene,
    quote = FALSE, sep = '\t'
  )

distinct(eqtlgen, SNP, SNPChr, SNPPos) %>%
  data.table::fwrite(
    .,
    file = file_eqtlgen_snp,
    quote = FALSE, sep = '\t'
  )

