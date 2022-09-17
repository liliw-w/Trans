############################################################
########### Check how Sigma from (1) independent null eQTLGen z-scores ###########
########### (2) DGN expression data ###########
########### affect PCO pvalues ###########
############################################################
library(data.table)
library(tidyverse)

module <- as.numeric(snakemake@params[['module']])
chr <- as.numeric(snakemake@params[['chr']])

thre_p_z <- 1e-4


### check how many snps will be used for this module under given null p threshold
file_null_SNP <- 'null_SNP/num_nullSNP.rds'
res_nullSNP <- readRDS(file_null_SNP)
res_nullSNP %>% filter(thre_z == thre_p_z & module == {{module}})



### READ DATA ###

#eqtlGen_transPCO_res_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/FDR/signals.qvalue.txt"
#eqtlGen_transPCO_res = fread(eqtlGen_transPCO_res_file, header = TRUE)

file_coexp_module <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds"
file_gene_meta <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/gene.meta.txt"
file.ex.var.regressed <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds"

file_eQTLGen_z <- paste0("sumStat_all/sumStat.all.module", module, ".rds")
file_indep_null_snp <- paste0("null_SNP/indep.SNP.module", module, ".", thre_p_z, ".txt")


### Filter out genes cross-map with the SNP's cis- genes
file_cross_map_genes <- "cross_map_genes.rds"
file_num_gene_snp_used_all<- 'p/num_gene_snp_used_module_all.Sigma_nullz.txt'


coexp_module <- readRDS(file_coexp_module)$moduleLabels; Nmodule <- max(coexp_module)
gene_meta <- fread(file_gene_meta, header = TRUE)
gene_meta$GeneNameConv <- sapply(gene_meta$GeneNameConv, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
datExpr <- readRDS(file.ex.var.regressed)

eQTLGen_z <- readRDS(file_eQTLGen_z)
indep_null_snp <- fread(file_indep_null_snp, header = FALSE)$V1
cross_map_genes <- readRDS(file_cross_map_genes)
num_gene_snp_used_all<- fread(file_num_gene_snp_used_all)


### Estimate the correlation matrix of zscore using extracted independent null z's, for the whole module ###
indep_null_z <- eQTLGen_z[eQTLGen_z$SNPmeta %in% indep_null_snp, -c(1:2)]
Sigma_null_z <- cor(indep_null_z)


### Apply PCO using this estimate Sigma ###
params1 <- "/project2/xuanyao/llw/DGN_no_filter_on_mappability/script/"
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))



p_all <- NULL
#num_gene_snp_used <- NULL
#for (chr in 1:22) {
  file_z <- paste0("/project2/xuanyao/llw/eQTLGen/z_dgn_166_module/z.module", module, ".chr", chr, ".txt.gz")
  z_mat <- fread(file_z)
  z_mat <- as.matrix(z_mat, rownames = TRUE)
  
  #if(ncol(z_mat)<2) next
  if(ncol(z_mat)<2){
  saveRDS(NULL, paste0('p/p.module', module, '.chr', chr, '.Sigma_nullz.rds'))
  q()
}

  ### Filter out genes on the same chr as the SNP
  gene_w_pos <- gene_meta[gene_meta$GeneNameConv %in% colnames(z_mat), ]
  gene_trans <- gene_w_pos[gene_w_pos$chr != paste0("chr", chr), ]$GeneNameConv
  
  z_mat_trans <- z_mat[, gene_trans, drop=FALSE]
  
  #if(ncol(z_mat_trans)<2) next
  if(ncol(z_mat_trans)<2){
  saveRDS(NULL, paste0('p/p.module', module, '.chr', chr, '.Sigma_nullz.rds'))
  q()
}

  
  ### run PCO for SNPs in two batches
  ## 1. SNPs share same Sigmas (gene used for PCO), vectorized PCO
  tmp_snp_share_sigma <- num_gene_snp_used_all %>%
    filter(module == {{module}} &
             SNPChr == {{chr}} &
             SNP %in% rownames(z_mat_trans) &
             n_gene_trans == n_gene_trans_cross) %>%
    pull(SNP)
  
  z_mat_trans_cross <- z_mat_trans
  ind_snp_shared_sigma <- rownames(z_mat_trans_cross) %in% tmp_snp_share_sigma
  
  ### Use Sigma from eQTLGen independent null zscores ###
  ind_sigma <- colnames(Sigma_null_z) %in% colnames(z_mat_trans_cross)
  Sigma <- Sigma_null_z[ind_sigma, ind_sigma, drop=FALSE]
  
  if(nrow(Sigma) >= 2){
    SigmaO <- ModifiedSigmaOEstimate(Sigma)
    p_snp <- ModifiedPCOMerged(Z.mat = z_mat_trans_cross[ind_snp_shared_sigma, , drop=FALSE],
                               Sigma = Sigma, SigmaO = SigmaO )
    p_all <- c(p_all, p_snp)
  }
  
  
  ## 2. SNPs with different genes, thus diff. signals, thus can only run PCO individually
  for(snp in rownames(z_mat_trans)[!ind_snp_shared_sigma] ){
    gene_cross <- cross_map_genes %>%
      filter(SNP %in% snp) %>%
      pull(cross_Geneid) %>%
      strsplit(split = ";") %>%
      unlist()
    
    z_mat_trans_cross <- z_mat_trans[, !colnames(z_mat_trans) %in% gene_cross, drop=FALSE]
    
    ### Use Sigma from eQTLGen independent null zscores ###
    ind_sigma <- colnames(Sigma_null_z) %in% colnames(z_mat_trans_cross)
    Sigma <- Sigma_null_z[ind_sigma, ind_sigma, drop=FALSE]
    
    if(nrow(Sigma) >= 2){
      SigmaO <- ModifiedSigmaOEstimate(Sigma)
      p_snp <- ModifiedPCOMerged(Z.mat = z_mat_trans_cross[rownames(z_mat_trans_cross) %in% snp, , drop=FALSE],
                                 Sigma = Sigma, SigmaO = SigmaO )
      p_all <- c(p_all, p_snp)
    }
    
    #num_gene_snp_used <- rbind(num_gene_snp_used,
    #                           c("module" = module,
    #                             "module_size" = sum(coexp_module == module),
    #                             "SNP" = snp,
    #                             "module_size_eqtlgen" = ncol(z_mat),
    #                             "n_gene_trans" = ncol(z_mat_trans),
    #                             "n_gene_trans_cross" = ncol(z_mat_trans_cross)) )
  }
  cat('Chr', chr, "has done!", "\n")
  
  
  res_p_all <- enframe(p_all, "SNP", "p") %>%
    left_join(cross_map_genes %>% select(1:4), by = "SNP")
  
  #res_num_gene_snp_used <- as_tibble(num_gene_snp_used) %>%
  #  left_join(cross_map_genes %>% select(1:4), by = "SNP") %>%
  #  mutate("num_nullSNP_indep" = res_nullSNP %>% filter(thre_z == thre_p_z & module == {{module}}) %>% pull(num_nullSNP_indep) )
  
  
  saveRDS(res_p_all, paste0('p/p.module', module, '.chr', chr, '.Sigma_nullz.rds'))
  
  
  #fwrite(res_num_gene_snp_used,
  #       paste0('p/num_gene_snp_used_module', module, '.Sigma_nullz.txt'),
  #       sep = "\t", quote = FALSE)
  
#}
cat('module', module, 'chr', chr, "has done!", "\n")




