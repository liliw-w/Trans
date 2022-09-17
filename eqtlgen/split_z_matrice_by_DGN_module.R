######################################################
######### split eQTLGen z matrices into (module, chr)-pair-based matrices #########
######### According to a known DGN module #########
######################################################
rm(list = ls())
library(data.table)
library(tidyverse)

file_dgn_module <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds'
file_dgn_gene_meta <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/gene.meta.txt'
file_eqtlGen_z <- list.files("/project2/xuanyao/llw/eQTLGen/", "^z.chr.*.txt.gz$")

dir_z <- 'z_dgn_166_module/'
file_num_snp_gene <- "z_dgn_166_module/num_snp_gene_z_module_chr.txt"


### read files
dgn_module <- readRDS(file_dgn_module)$moduleLabels
dgn_gene_meta <- fread(file_dgn_gene_meta, header = TRUE)
tmp_z <- fread(file_eqtlGen_z[3], header = TRUE)


### what are the DGN genes and eQTLGen genes (take one chr for example here)
dgn_gene <- dgn_gene_meta %>%
  filter(gene %in% names(dgn_module)) %>%
  separate(GeneNameConv, c("gene_ENSG", NA), sep = "[.]") %>%
  select(1:5) %>%
  mutate("module" = dgn_module[gene]) %>%
  arrange(module)

Nmodule <- max(dgn_module)

eqtlGen_gene <- colnames(tmp_z)


### plot the DGN genes and their overlap with eQTLGen genes in each modules
ggplot(data = dgn_gene %>% mutate("if_eqtlGen" = gene_ENSG %in% eqtlGen_gene),
       aes(x = factor(module), fill = if_eqtlGen)) +
  geom_bar()


### split z matrices
num_snp_gene = tibble()
for(chr in 1:22){
  z_chr_mat <- fread(paste0("z.chr", chr, ".txt.gz"), header = TRUE)
  z_chr_mat <- z_chr_mat %>% filter( complete.cases(.) )
  
  for(m in 1:Nmodule){
    tmp_dgn_gene <- filter(dgn_gene, module == m)
    z_col_name <- c('snp', tmp_dgn_gene$gene_ENSG )
    z_module_chr_mat <- z_chr_mat %>% select(any_of( z_col_name ))
    
    
    fwrite(z_module_chr_mat,
           paste0(dir_z, "z.module", m, ".chr", chr, ".txt.gz"),
           quote = FALSE, sep = "\t")
    cat("module:", m, "\n")
    
    num_snp_gene = rbind(num_snp_gene,
                             c(m, chr, nrow(z_module_chr_mat), ncol(z_module_chr_mat)-1 ))
  }
  
  cat(chr, "\n")
}

### keep track of how many eQTLGen SNP and genes are left in each z matrices
colnames(num_snp_gene) <- c("module", "chr", "num_of_snp", "num_of_gene")
fwrite(num_snp_gene,
       file_num_snp_gene,
       quote = FALSE, sep = "\t")

