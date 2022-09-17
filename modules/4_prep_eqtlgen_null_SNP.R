############################################################
####### Extract z scores from eQTLGen sum stat #######
####### and assemble them into pre-defined gene modules #######
####### and pick out independent SNPs #######
####### for estimating \Sigma #######
############################################################
# load packages -----
rm(list = ls())
library(data.table)


# I/O & paras -----
file_coexp_module <- '/project2/xuanyao/llw/eQTLGen/MSigDB/coexp.module.rds'
file_gene_meta <- '/project2/xuanyao/llw/eQTLGen/eqtlgen_genes.txt'
file_eqtlgen <- '/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz'
script_indep_SNP <- '/home/liliw1/Trans/eqtlgen/indep.SNP.sh'

dir_trans_sumstat <- "sumStat_trans/"
dir_all_sumstat <- "sumStat_all/"
dir_null_snp <- "null_SNP/"

## output -----
dir.create(dir_trans_sumstat, showWarnings = FALSE)
dir.create(dir_all_sumstat, showWarnings = FALSE)
dir.create(dir_null_snp, showWarnings = FALSE)

file_res_nullSNP <- "null_SNP/num_nullSNP.rds"


# read files -----
gene_meta <- fread(file_gene_meta)
coexp_module <- readRDS(file_coexp_module)$moduleLabels
eqtlgen <- fread(file_eqtlgen,
                 select = c('SNP', 'SNPChr', 'SNPPos', 'Zscore', 'Gene', 'GeneSymbol'))


# organize data -----
## convert the names of genes in modules to ENSG id -----
names(coexp_module) <- gene_meta[match(names(coexp_module), gene_meta$GeneSymbol), Gene]
Nmodule <- max(coexp_module)

## remove dots in ENSG gene names -----
names(coexp_module) <- sapply(names(coexp_module), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])


# 1. Extract (trans) summary stat for each modules from the eQTLGen sum stat -----
for(module in 1:Nmodule){
  gene_module <- names(coexp_module)[coexp_module==module]
  eqtlgen_module <- eqtlgen[eqtlgen$Gene %in% gene_module, ]
  
  saveRDS(
    eqtlgen_module,
    paste0(dir_trans_sumstat, "sumStat.trans.module", module, ".rds")
  )
  cat("module:", module, "\n")
}


# 2. Combine trans and cis summary stat for each module -----
## And extract the common SNPs, i.e. SNPs with summary stat for all genes within each module.
for(module in 1:Nmodule){
  eqtlgen_trans_module <- readRDS(paste0(dir_trans_sumstat, "sumStat.trans.module", module, ".rds"))
  
  eqtlgen_trans_module$SNPmeta <- paste(eqtlgen_trans_module$SNPChr, eqtlgen_trans_module$SNPPos, sep = ":")
  
  ab_z <- dcast(eqtlgen_trans_module, SNP+SNPmeta~Gene, fun.aggregate = max,
                value.var = 'Zscore',
                drop = TRUE, fill = NA)
  ab_z_na <- is.na(ab_z)
  rem_ind <- rowSums(ab_z_na[, -c(1:2)])
  ab_z_remain <- ab_z[rem_ind == 0, ]
  
  
  saveRDS(
    ab_z_remain,
    paste0(dir_all_sumstat, "sumStat.all.module", module, ".rds")
  )
  
  cat("module:", module, "\n")
}


# 3. Indep null SNPs (and #Uniq null SNPs) for each module under various z-score threshold for being null -----
res_nullSNP <- NULL
for(module in 1:Nmodule){
  ab_z_remain <- readRDS(paste0(dir_all_sumstat, "sumStat.all.module", module, ".rds"))
  NSNP <- nrow(ab_z_remain); module_size <- ncol(ab_z_remain)-2
  
  for(thre_p_z in c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)){
    thre_z2 <- qchisq(1-thre_p_z, 1)
    altSNP_ind <- apply((ab_z_remain[, -c(1:2)])^2 > thre_z2, 1, sum) > 0
    
    file_uniq <- paste0(dir_null_snp, "unique.SNP.module", module, ".", thre_p_z, ".txt")
    file_indep <- paste0(dir_null_snp, "indep.SNP.module", module, ".", thre_p_z, ".txt")
    fwrite(ab_z_remain[!altSNP_ind, "SNPmeta"], file_uniq,
           quote = FALSE, sep = "\t", col.names = FALSE)
    command <- paste("bash", script_indep_SNP, file_uniq, file_indep)
    system(command)
    num_nullSNP_indep <- as.numeric(system(paste0("cat ", file_indep, " | wc -l"), intern = TRUE))
    
    res_nullSNP <- rbind(res_nullSNP, data.table("thre_z" = thre_p_z,
                                                 "num_nullSNP_uniq" = NSNP - sum(altSNP_ind),
                                                 "num_nullSNP_indep" = num_nullSNP_indep,
                                                 "module" = module,
                                                 "module_size" = module_size))
  }
}


# print out key message and write out -----
saveRDS(res_nullSNP, file_res_nullSNP)

