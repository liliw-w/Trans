##############################################
########### extract z for genes in module m and corresponding snps ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    /project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt
    /project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds
    /project2/xuanyao/llw/DGN_no_filter_on_mappability/result/gene.meta.txt
    4
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

files_signal <- args[1]
file_gene_module <- args[2]
file_gene_meta <- args[3]
m <- as.numeric(args[4])


## output -----


# read files -----
signal <- data.table::fread(files_signal, col.names = c("sig", "p", "q"))
gene_module <- readRDS(file_gene_module)$moduleLabels %>%
  enframe(name = 'gene', value = 'module')
gene_meta <- data.table::fread(file_gene_meta, header = TRUE)


# extract signals and module genes for module m -----
signal <- separate(
  signal, 
  col = 'sig', 
  into = c('module', 'chr_snp', 'chr_pos'), 
  convert = TRUE
) %>%
  filter(module == str_glue('module{m}'))

gene_module <- filter(gene_module, module == !!m)


# extract z for module m and corresponding snps -----
signal_z <- lapply(unique(signal$chr_snp), function(chr) {
  signal_chr = filter(signal, chr_snp == !!chr)
  
  file_z = list.files(
    path = '/project2/xuanyao/llw/DGN_no_filter_on_mappability/z',
    pattern = str_glue('z.module{m}.chr{chr}.txt.gz'),
    full.names = TRUE
  )
  
  z = data.table::fread(file_z)
  
  cat(str_glue('{chr} \n'))
  
  return(
    filter(z, snp %in% str_glue('{signal_chr$chr_snp}:{signal_chr$chr_pos}'))
  )
})


out <- bind_rows(signal_z) %>%
  # make genes as rows, snps as cols -----
  pivot_longer(cols = !snp, names_to = "gene", values_to = 'z') %>%
  pivot_wider(names_from = 'snp', values_from = 'z') %>%
  
  # add gene meta -----
  left_join(
    select(gene_meta, gene, chr),
    by = 'gene'
  ) %>%
  relocate(
    chr, .after = 'gene'
  ) %>%
  rename('gene_chr' = 'chr')


# save files -----
saveRDS(signal, str_glue('signal_module{m}.rds'))
saveRDS(tmp1, str_glue('z_module{m}.rds'))

