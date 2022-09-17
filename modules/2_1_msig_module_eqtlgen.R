##############################################
########### Re-format MSigDB modules using eQTLGen genes ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')


# I/O & paras -----
file_msigdb_module_df <- '/project2/xuanyao/llw/MODULES/MSigDB/h.all.v7.4.symbols.txt'
file_eqtlgen_genes <- "/project2/xuanyao/llw/eQTLGen/eqtlgen_genes.txt"

## output -----
file_coexp_module <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/coexp.module.rds'
file_msigdb_module_eqtlgen <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/msigdb_module_eqtlgen.txt'


# read files -----
msigdb_module <- fread(file_msigdb_module_df, header = TRUE)
eqtlgen_genes <- fread(file_eqtlgen_genes)


# organize data -----
## overlap info with eqtlgen genes and module size -----
msigdb_module <- msigdb_module %>%
  mutate(if_eqtlgen = gene %in% !!eqtlgen_genes$GeneSymbol,
         anno = category,
         category = fct_infreq(category)) %>%
  group_by(anno) %>%
  mutate(
    module_size = n(),
    module_size_eqtlgen = sum(if_eqtlgen),
  ) %>%
  ungroup()


## recode module annot to numbers -----
msigdb_module$category <- factor(
  msigdb_module$category,
  labels = 1:nlevels(msigdb_module$category)
)

## order rows by module id and if in eqtlgen -----
msigdb_module <- arrange(msigdb_module, category, desc(if_eqtlgen))


## format modules as input for the pipeline -----
coexp_module <- as.numeric(msigdb_module$category) %>%
  setNames(msigdb_module$gene) %>%
  list("moduleLabels" = .)

coexp_module$if_eqtlgen <- msigdb_module$if_eqtlgen
names(coexp_module$if_eqtlgen) <- msigdb_module$gene

coexp_module$moduleName <- select(msigdb_module, anno, category) %>%
  distinct() %>%
  mutate(category = as.numeric(category)) %>%
  deframe()


# print out key message and write out -----
saveRDS(
  coexp_module,
  file_coexp_module
)
fwrite(
  msigdb_module,
  file_msigdb_module_eqtlgen,
  quote = FALSE, sep = "\t"
)

cat(
  "There are", nrow(eqtlgen_genes), "genes analyzed in eQTLGen for trans- associations. \n\n",
  
  "There are", nrow(msigdb_module), "(uniquely", distinct(msigdb_module, gene) %>% nrow(), ") genes in MSigDB modules. \n\n",
  
  "And the overlapped genes are", sum(msigdb_module$gene %in% eqtlgen_genes$GeneSymbol),
  "(uniquely", sum(unique(msigdb_module$gene) %in% eqtlgen_genes$GeneSymbol), ") genes. \n\n"
)
