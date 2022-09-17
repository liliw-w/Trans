############################################################
########### Check for each SNP, how many genes are remained in each module ###########
########### after remove cross-map genes of their cis-genes ###########
########### i.e. final genes used for PCO ###########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)


# I/O & paras -----
file_gene_meta <- '/project2/xuanyao/llw/eQTLGen/eqtlgen_genes.txt'
file_cross_map_genes <- "/project2/xuanyao/llw/eQTLGen_est_Sigma/cross_map_genes.rds"
file_coexp_module <- '/project2/xuanyao/llw/eQTLGen/MSigDB/coexp.module.rds'
file_msigdb_module_eqtlgen <- '/project2/xuanyao/llw/eQTLGen/MSigDB/msigdb_module_eqtlgen.txt'
file_null_SNP <- 'null_SNP/num_nullSNP.rds'

thre_p_z <- 1e-4

## output -----
file_num_gene_snp_used <- 'p/num_gene_snp_used_module_all.Sigma_nullz.txt'


# read files -----
gene_meta <- fread(file_gene_meta, header = TRUE)
# to remove cross-map genes with cis-genes of each SNP
cross_map_genes <- readRDS(file_cross_map_genes)
coexp_module <- readRDS(file_coexp_module)$moduleLabels; Nmodule <- max(coexp_module)
msigdb_module_eqtlgen <- fread(file_msigdb_module_eqtlgen, header = TRUE)
# check how many snps will be used for this module under given null p threshold
res_nullSNP <- readRDS(file_null_SNP)


# organize data -----
## remove dot in ENSG gene names -----
gene_meta$GeneNameConv <- sapply(gene_meta$Gene, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])



# for (module, chr) pair, how many genes will be used for each SNP, after removing genes on same chr and cross map genes with cis genes -----
num_gene_snp_used <- NULL
for(module in 1:Nmodule){
  module_size_eqtlgen <- distinct(msigdb_module_eqtlgen, category, module_size_eqtlgen) %>%
    filter(category == !!module) %>%
    pull(module_size_eqtlgen)
  
  for (chr in 1:22) {
    file_z <- paste0("/project2/xuanyao/llw/eQTLGen/MSigDB/z/z.module", module, ".chr", chr, ".txt.gz")
    z_mat <- fread(file_z)
    z_mat <- as.matrix(z_mat, rownames = TRUE)
    
    ## Filter out genes in the module on the same chr as the SNP -----
    gene_w_pos <- gene_meta[gene_meta$GeneNameConv %in% colnames(z_mat), ]
    gene_trans <- gene_w_pos[gene_w_pos$GeneChr != paste0("chr", chr), ]$GeneNameConv
    
    z_mat_trans <- z_mat[, gene_trans, drop=FALSE]
    
    
    ## Filter out genes in the module cross map with cis genes of the snp -----
    for(snp in rownames(z_mat_trans) ){
      gene_cross <- cross_map_genes %>%
        filter(SNP %in% snp) %>%
        pull(cross_Geneid) %>%
        strsplit(split = ";") %>%
        unlist()
      
      z_mat_trans_cross <- z_mat_trans[, !colnames(z_mat_trans) %in% gene_cross, drop=FALSE]
      
      num_gene_snp_used <- rbind(num_gene_snp_used,
                                 c("module" = module,
                                   "module_size" = sum(coexp_module == module),
                                   "SNP" = snp,
                                   "module_size_eqtlgen" = module_size_eqtlgen,
                                   "n_gene_trans" = ncol(z_mat_trans),
                                   "n_gene_trans_cross" = ncol(z_mat_trans_cross)) )
    }
    cat('Chr', chr, "has done!", "\n")
  }
  cat('module', module, "has done!", "\n")
  
}

## add SNP meta info -----
num_gene_snp_used <- as_tibble(num_gene_snp_used) %>%
  left_join(cross_map_genes %>% select(1:4), by = "SNP") %>%
  type_convert(col_types = cols(meta = col_character())) %>%
  left_join(
    filter(res_nullSNP, thre_z == thre_p_z) %>% select(module, num_nullSNP_indep),
    by = "module"
  )


# print out key message and write out -----
fwrite(
  num_gene_snp_used,
  file_num_gene_snp_used,
  sep = "\t",
  quote = FALSE
)



