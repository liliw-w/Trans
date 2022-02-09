rm(list = ls())
library(data.table)
library(tidyverse)

file.qtl = 'FDR/signals.chr.module.perm10.fdr10.txt'
file.gene.meta = '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
file_signal_cis_genes = 'postanalysis/signal_cis_genes.txt'

qtl = fread(file.qtl, header = FALSE, col.names = c("signal", "p", "q"))
gene.meta = fread(file.gene.meta, header = TRUE)

qtl = qtl %>%
  separate("signal",
           into = c("module", "chr", "pos"),
           sep = ":", remove = FALSE, convert = TRUE) %>%
  unite(col = "SNP",
        c("chr", "pos"),
        sep = ":", remove = FALSE)

gene.meta = gene.meta %>%
  filter(Class %in% c("protein_coding", "lincRNA") & Chromosome %in% paste0("chr",1:22)) %>%
  separate("Chromosome", c(NA, "chr"), sep = "chr", remove = FALSE)

cis.gene.meta = sapply(1:nrow(qtl), function(x){
  tmp_qtl = qtl[x, ]
  tmp_cis.gene.meta = gene.meta %>%
    filter(chr %in% tmp_qtl$chr) %>%
    mutate("dis" = abs(tmp_qtl$pos - Start)) %>%
    filter(dis < 5e+4) %>%
    arrange(dis)
  c(
    tmp_cis.gene.meta$GeneSymbol[1],
    tmp_cis.gene.meta$dis[1],
    paste(tmp_cis.gene.meta$GeneSymbol, collapse = ";"),
    paste(tmp_cis.gene.meta$dis, collapse = ";")
  )
})
cis.gene.meta = as.data.table(t(cis.gene.meta))
colnames(cis.gene.meta) = c("nearest_gene", "nearest_dis", "near_genes", "near_dis")


qtl = bind_cols(qtl, cis.gene.meta)

fwrite(qtl, file_signal_cis_genes, quote = FALSE, sep = "\t")
