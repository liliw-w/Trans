##############################################
########### eQTLGen specific signals not in DGN ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)



# paras and I/O -----
file_eqtlgen_sig_dgn <- 'postanalysis/eqtlgen_in_dgn.txt'
file_dgn_sig <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'

file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
dis_cis <- 1e+6

buffer_reg <- c(1e+5, 2e+5, 5e+5, 1e+6)

file_out <- 'postanalysis/eqtlgen_spec_to_dgn.txt'


# read files -----
eqtlgen_sig_dgn <- fread(file_eqtlgen_sig_dgn, header = TRUE)
dgn_sig <- fread(file_dgn_sig, header = FALSE, col.names = c("signal", "p", "q") )
gene_meta <- fread(file_gene_meta, header = TRUE)


# organize data -----
dgn_sig <- dgn_sig %>%
  separate(signal,
           into = c('module', 'chr', 'pos'),
           sep = ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_id", c(chr, pos), sep = ":", remove = FALSE)


# eQTLGen specific signals -----
# defined as SNPs not within cis- region of any DGN signals
min_dis_to_dgn_sig <- sapply(1:nrow(eqtlgen_sig_dgn), function(x){
  tmp_df = eqtlgen_sig_dgn[x, ]
  
  dgn_sig %>%
    filter(module == tmp_df$module & chr == tmp_df$SNPChr) %>%
    mutate("dis" = abs(pos - tmp_df$SNPPos)) %>%
    pull(dis) %>%
    min()
})

# add distance of each eqtlgen sig to nearest dgn sig
eqtlgen_sig_dgn <- mutate(eqtlgen_sig_dgn, "min_dis_to_dgn_sig" = min_dis_to_dgn_sig)

# add specificity of each row
eqtlgen_sig_dgn <- cbind(
  eqtlgen_sig_dgn,
  
  map2_dfc(
    eqtlgen_sig_dgn %>% select(min_dis_to_dgn_sig),
    buffer_reg,
    ~.x<=.y
  ) %>%
    rename_with(~paste0("if_spec_", buffer_reg))
)


# add cis- genes -----
# use only "protein_coding" & "lincRNA" as cis- genes
gene_meta <- gene_meta %>%
  filter(Class %in% c("protein_coding", "lincRNA") & Chromosome %in% paste0("chr",1:22)) %>%
  separate("Chromosome", c(NA, "chr"), sep = "chr", remove = FALSE, convert = TRUE)

# cis- genes of every snp
cis_gene_meta <- sapply(1:nrow(eqtlgen_sig_dgn), function(x){
  tmp_qtl = eqtlgen_sig_dgn[x, ]
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

# bind cis- gene columns to signals df
eqtlgen_sig_dgn <- bind_cols(eqtlgen_sig_dgn, cis_gene_meta)


# save eQTLGen unique signals -----
fwrite(
  eqtlgen_sig_dgn,
  file_out,
  quote = FALSE, sep = "\t"
)


cat(
  "Out of", nrow(eqtlgen_sig_dgn), "eQTLGen (SNP, module) signal pairs that were also analyze by DGN, \n\n",
  
  colSums(select(eqtlgen_sig_dgn, starts_with("if_spec"))), "pairs are near at least one DGN signal within same loci of distance: \n\n",
  
  buffer_reg, "base pairs, respectively. \n\n",
  
  sum(eqtlgen_sig_dgn$min_dis_to_dgn_sig == Inf), "pairs are not close to any DGN signals. \n\n"
)

