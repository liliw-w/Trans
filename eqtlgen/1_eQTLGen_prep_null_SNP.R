############################################################
####### Extract z scores from eQTLGen sum stat #######
####### and assemble them into DGN gene modules #######
####### and pick out independent SNPs #######
####### for estimating \Sigma #######
############################################################

## load packages
rm(list = ls())
library(data.table)

## I/O
file_meta <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/gene.meta.txt'
file_module_info <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds'
file_eqtlgen <- '/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz'
script_indep_SNP <- '/project2/xuanyao/llw/eQTLGen_est_Sigma/indep.SNP.sh'


## read gene meta file and module file
meta <- fread(file_meta)
module_info <- readRDS(file_module_info)$moduleLabels

## convert the names of genes in modules to ENSG id
names(module_info) <- meta[match(names(module_info), meta$gene), GeneNameConv]
Nmodule <- max(module_info)

## remove dots in ENSG gene names
names(module_info) <- sapply(names(module_info), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])


## 1. Extract (trans) summary stat for each modules from the eQTLGen sum stat
eqtlgen <- fread(file_eqtlgen,
          select = c('SNP', 'SNPChr', 'SNPPos', 'Zscore', 'Gene', 'GeneSymbol'))
for(module in 1:Nmodule){
  gene_module <- names(module_info)[module_info==module]
  eqtlgen_module <- eqtlgen[eqtlgen$Gene %in% gene_module, ]

  saveRDS(eqtlgen_module, paste0("sumStat_trans/sumStat.trans.module", module, ".rds"))
  cat("module:", module, "\n")
}


## Extract (cis) summary stat for each modules
#a = fread('/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz',
#          select = c('SNP', 'SNPChr', 'SNPPos', 'Zscore', 'Gene', 'GeneSymbol'))
#for(module in 1:Nmodule){
#  module1 = names(module_info)[module_info==module]
#  a_module1 = a[a$Gene %in% module1, ]
#
#  saveRDS(a_module1, paste0("sumStat_cis/sumStat.cis.module", module, ".rds"))
#  cat("module:", module, "\n")
#}


## 2. Combine trans and cis summary stat for each module
## And extract the common SNPs, i.e. SNPs with summary stat for all genes within each module.
for(module in 1:Nmodule){
  eqtlgen_cis_module <- NULL #readRDS(paste0("sumStat_cis/sumStat.cis.module", module, ".rds"))
  eqtlgen_trans_module <- readRDS(paste0("sumStat_trans/sumStat.trans.module", module, ".rds"))
  
  eqtlgen_both_module <- rbind(eqtlgen_cis_module, eqtlgen_trans_module)
  eqtlgen_both_module$SNPmeta <- paste(eqtlgen_both_module$SNPChr, eqtlgen_both_module$SNPPos, sep = ":")

  ab_z <- dcast(eqtlgen_both_module, SNP+SNPmeta~Gene, fun.aggregate = max,
               value.var = 'Zscore',
               drop = TRUE, fill = NA)
  ab_z_na <- is.na(ab_z)
  rem_ind <- rowSums(ab_z_na[, -c(1:2)])
  ab_z_remain <- ab_z[rem_ind == 0, ]
  
  
  saveRDS(ab_z_remain, paste0("sumStat_all/sumStat.all.module", module, ".rds"))
  
  cat("module:", module, "\n")
}


## 3. Indep null SNPs (and #Uniq null SNPs) for each module under various z-score threshold for being null
res_nullSNP <- NULL
for(module in 1:Nmodule){
  ab_z_remain <- readRDS(paste0("sumStat_all/sumStat.all.module", module, ".rds"))
  NSNP <- nrow(ab_z_remain); module_size <- ncol(ab_z_remain)-2

  for(thre_p_z in c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)){
    thre_z2 <- qchisq(1-thre_p_z, 1)
    altSNP_ind <- apply((ab_z_remain[, -c(1:2)])^2 > thre_z2, 1, sum) > 0

    file_uniq <- paste0("null_SNP/unique.SNP.module", module, ".", thre_p_z, ".txt")
    file_indep <- paste0("null_SNP/indep.SNP.module", module, ".", thre_p_z, ".txt")
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
saveRDS(res_nullSNP, "null_SNP/num_nullSNP.rds")


