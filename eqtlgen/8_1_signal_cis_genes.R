##############################################
########### nearest and cis genes of signals ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
ratio <- 50

file_qtl <- paste0('postanalysis/signal_rm_infl_ratio_', ratio, '.txt')
file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
dis_cis <- 1e+6

## output -----
file_signal_cis_genes <- paste0('postanalysis/signal_cis_genes_rm_infl_ratio_', ratio, '.txt')


# read files -----
qtl <- fread(file_qtl, header = TRUE)
gene_meta <- fread(file_gene_meta, header = TRUE)


# organize data -----
## extract signals' module, chr, pos -----
#qtl <- qtl %>%
#  separate("signal",
#           into = c("module", "chr", "pos"),
#           sep = ":", remove = FALSE, convert = TRUE) %>%
#  unite(col = "SNP",
#        c("chr", "pos"),
#        sep = ":", remove = FALSE)

## extract protein_coding, lincRNA, auto-chr genes -----
gene_meta <- gene_meta %>%
  filter(Class %in% c("protein_coding", "lincRNA") & Chromosome %in% paste0("chr",1:22)) %>%
  separate("Chromosome", c(NA, "chr"), sep = "chr", remove = FALSE)


# nearest and cis genes of each signal, with distance -----
cis_gene_meta <- sapply(1:nrow(qtl), function(x){
  tmp_qtl = qtl[x, ]
  tmp_cis_gene_meta = gene_meta %>%
    filter(chr %in% tmp_qtl$SNPChr) %>%
    mutate("dis" = abs(tmp_qtl$SNPPos - Start)) %>%
    filter(dis < dis_cis/2) %>%
    arrange(dis)
  c(
    tmp_cis_gene_meta$GeneSymbol[1],
    tmp_cis_gene_meta$dis[1],
    paste(tmp_cis_gene_meta$GeneSymbol, collapse = ";"),
    paste(tmp_cis_gene_meta$dis, collapse = ";")
  )
})
cis_gene_meta <- as.data.table(t(cis_gene_meta))
colnames(cis_gene_meta) <- c("nearest_gene", "nearest_dis", "near_genes", "near_dis")


## add cis gene info to signals -----
qtl <- bind_cols(qtl, cis_gene_meta)



# print out key message or write out -----
fwrite(
  qtl,
  file_signal_cis_genes,
  quote = FALSE,
  sep = "\t"
)

