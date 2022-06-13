##############################################
########### cis- genes of eQTLGen trans signals ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)



# paras and I/O -----
dis_cis1 <- 1e+5
dis_cis2 <- 1e+6

file_eqtlgen_sig <- list.files('postanalysis', '^signal_rm_infl_ratio_\\d+[.]txt$', full.names = TRUE)
file_gene_anno <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'

file_signal_cis_genes <- paste0('postanalysis/', basename(file_eqtlgen_sig), '_cis_genes.txt')


# read files -----
eqtlgen_sig <- fread(file_eqtlgen_sig, header = TRUE)
gene_anno <- fread(file_gene_anno, header = TRUE)


# organize data -----
## keep only protein coding or lincRNA genes
gene_anno <- gene_anno %>%
  filter(
    Class %in% c("protein_coding", "lincRNA") &
      Chromosome %in% paste0("chr",1:22)
  ) %>%
  separate(
    Geneid,
    c('Geneid', NA),
    sep = "[.]"
  ) %>%
  separate(
    "Chromosome",
    c(NA, "chr"),
    sep = "chr",
    remove = FALSE,
    convert = TRUE
  )



# cis genes of each snp -----
## record three types of near genes: nearest, genes within dis1, genes within dis2
cis_gene <- lapply(1:nrow(eqtlgen_sig), function(x){
  tmp_qtl = eqtlgen_sig[x, ]
  
  tmp_cis_gene = gene_anno %>%
    filter(chr == tmp_qtl$SNPChr) %>%
    mutate("dis" = abs(tmp_qtl$SNPPos - Start)) %>%
    filter(dis < max(dis_cis1, dis_cis2)/2) %>%
    arrange(dis) %>%
    mutate(
      "if_in_dis1" = dis < dis_cis1/2,
      "if_in_dis2" = dis < dis_cis2/2
    )
  
  c(
    tmp_cis_gene$GeneSymbol[1],
    tmp_cis_gene$Geneid[1],
    tmp_cis_gene$dis[1],
    paste(tmp_cis_gene[tmp_cis_gene$if_in_dis1, GeneSymbol], collapse = ";"),
    paste(tmp_cis_gene[tmp_cis_gene$if_in_dis1, Geneid], collapse = ";"),
    paste(tmp_cis_gene[tmp_cis_gene$if_in_dis1, dis], collapse = ";"),
    paste(tmp_cis_gene[tmp_cis_gene$if_in_dis2, GeneSymbol], collapse = ";"),
    paste(tmp_cis_gene[tmp_cis_gene$if_in_dis2, Geneid], collapse = ";"),
    paste(tmp_cis_gene[tmp_cis_gene$if_in_dis2, dis], collapse = ";")
  )
}
)

## re-format and add col names
cis_gene <- rbindlist(list(cis_gene)) %>% t() %>% as_tibble()
colnames(cis_gene) <- c(
  "nearest_gene", "nearest_gene_symbol", "nearest_dis",
  paste0("near_genes_", dis_cis1),
  paste0("near_genes_symbol_", dis_cis1),
  paste0("near_dis_", dis_cis1),
  paste0("near_genes_", dis_cis2),
  paste0("near_genes_symbol_", dis_cis2),
  paste0("near_dis_", dis_cis2)
)


## combine the cis genes to the original signal df
eqtlgen_sig_cis <- bind_cols(eqtlgen_sig, cis_gene) %>%
  arrange(SNPChr, SNPPos, p)




# print out key message and write out  -----
fwrite(eqtlgen_sig_cis, file_signal_cis_genes,
       quote = FALSE, sep = "\t")

