##############################################
########### split eQTLGen z matrices into (module, chr)-pair-based matrices ###########
########### According to MSigDB modules ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_coexp_module <- '/project2/xuanyao/llw/eQTLGen/MSigDB/coexp.module.rds'
file_gene_meta <- '/project2/xuanyao/llw/eQTLGen/eqtlgen_genes.txt'
file_eqtlGen_z <- list.files(
  "/project2/xuanyao/llw/eQTLGen",
  "^z.chr\\d+.txt.gz$",
  full.names = TRUE
)

## output -----
dir_z <- '/project2/xuanyao/llw/eQTLGen/MSigDB/z/'
file_num_snp_gene <- "/project2/xuanyao/llw/eQTLGen/MSigDB/num_snp_gene_z_module_chr.txt"


# read files -----
coexp_module <- readRDS(file_coexp_module)$moduleLabels
gene_meta <- fread(file_gene_meta, header = TRUE)


# organize data -----
## genes in modules -----
coexp_module_gene <- enframe(coexp_module, name = "GeneSymbol", value = "module") %>%
  inner_join(distinct(gene_meta, GeneSymbol, .keep_all = TRUE), by = "GeneSymbol") %>%
  arrange(module)


## number of modules -----
Nmodule <- max(coexp_module)



# split z matrices -----
num_snp_gene <- tibble()
for(file_eqtlGen_z_chr in file_eqtlGen_z){
  chr <- str_extract(basename(file_eqtlGen_z_chr), "\\d+")
  
  z_chr_mat <- fread(file_eqtlGen_z_chr, header = TRUE)
  z_chr_mat <- z_chr_mat %>% filter( complete.cases(.) )
  
  
  for(m in 1:Nmodule){
    tmp_module_gene <- filter(coexp_module_gene, module == m)
    z_col_name <- c('snp', tmp_module_gene$Gene )
    z_module_chr_mat <- z_chr_mat %>% select(any_of( z_col_name ))
    
    
    fwrite(z_module_chr_mat,
           paste0(dir_z, "z.module", m, ".chr", chr, ".txt.gz"),
           quote = FALSE, sep = "\t")
    
    cat("module:", m, "\n")
    
    
    num_snp_gene <- rbind(num_snp_gene, 
                          c(m, chr, nrow(z_module_chr_mat), ncol(z_module_chr_mat)-1 ))
  }
  
  cat("chr:", chr, "\n")
}


# keep track of how many eQTLGen SNP and genes are left in each z matrices -----
colnames(num_snp_gene) <- c("module", "chr", "num_of_snp", "num_of_gene")
fwrite(num_snp_gene,
       file_num_snp_gene,
       quote = FALSE, sep = "\t")

