############################################################
########### Check for each SNP, how many genes are remained in each module ###########
########### after remove cross-map genes of their cis-genes ###########
########### i.e. final genes used for PCO ###########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)

thre_p_z <- 1e-4


### READ DATA ###
file_null_SNP <- 'null_SNP/num_nullSNP.rds'
file_coexp_module <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds"
file_gene_meta <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/gene.meta.txt"
file_cross_map_genes <- "cross_map_genes.rds"


# check how many snps will be used for this module under given null p threshold
res_nullSNP <- readRDS(file_null_SNP)
coexp_module <- readRDS(file_coexp_module)$moduleLabels; Nmodule <- max(coexp_module)
gene_meta <- fread(file_gene_meta, header = TRUE)
gene_meta$GeneNameConv <- sapply(gene_meta$GeneNameConv, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
# to remove cross-map genes with cis-genes of each SNP
cross_map_genes <- readRDS(file_cross_map_genes)


num_gene_snp_used <- NULL
for(module in 1:Nmodule){
  for (chr in 1:22) {
    file_z <- paste0("/project2/xuanyao/llw/eQTLGen/z_dgn_166_module/z.module", module, ".chr", chr, ".txt.gz")
    z_mat <- fread(file_z)
    z_mat <- as.matrix(z_mat, rownames = TRUE)
    
    ### Filter out genes on the same chr as the SNP
    gene_w_pos <- gene_meta[gene_meta$GeneNameConv %in% colnames(z_mat), ]
    gene_trans <- gene_w_pos[gene_w_pos$chr != paste0("chr", chr), ]$GeneNameConv
    
    z_mat_trans <- z_mat[, gene_trans]
    
    
    for(snp in rownames(z_mat_trans) ){
      gene_cross <- cross_map_genes %>%
        filter(SNP %in% snp) %>%
        pull(cross_Geneid) %>%
        strsplit(split = ";") %>%
        unlist()
      
      z_mat_trans_cross <- z_mat_trans[, !colnames(z_mat_trans) %in% gene_cross]
      
      num_gene_snp_used <- rbind(num_gene_snp_used,
                                 c("module" = module,
                                   "module_size" = sum(coexp_module == module),
                                   "SNP" = snp,
                                   "module_size_eqtlgen" = ncol(z_mat),
                                   "n_gene_trans" = ncol(z_mat_trans),
                                   "n_gene_trans_cross" = ncol(z_mat_trans_cross)) )
    }
    cat('Chr', chr, "has done!", "\n")
  }
  cat('module', module, "has done!", "\n")
  
}

### add SNP meta info
res_num_gene_snp_used <- as_tibble(num_gene_snp_used) %>%
  left_join(cross_map_genes %>% select(1:4), by = "SNP") %>%
  mutate("num_nullSNP_indep" = res_nullSNP %>% filter(thre_z == thre_p_z & module == {{module}}) %>% pull(num_nullSNP_indep) )

### save result
fwrite(res_num_gene_snp_used,
       'p/num_gene_snp_used_module_all.Sigma_nullz.txt',
       sep = "\t", quote = FALSE)



