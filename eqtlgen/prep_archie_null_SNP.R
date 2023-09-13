############################################################
####### Extract z scores from eQTLGen sum stat #######
####### and assemble them into gene modules #######
####### and pick out independent SNPs #######
####### for estimating \Sigma #######
############################################################
# load packages -----
rm(list = ls())
library(tidyverse)


# I/O & paras -----
file_signal_archie <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/archie/archie_results.xlsx'
file_eqtlgen <- '/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz'
script_indep_SNP <- '~/Trans/eqtlgen/indep.SNP.sh'


## output -----
file_module_archie <- 'null_SNP/module_archie.rds'
file_z <- "null_SNP/sumStat_archie.rds"
file_num_nullsnp <- "null_SNP/num_nullSNP.rds"


# read files -----
eqtlgen <- data.table::fread(
  file_eqtlgen,
  select = c('SNP', 'SNPChr', 'SNPPos', 'Zscore', 'Gene', 'GeneSymbol')
)
signal_archie <- readxl::read_excel(file_signal_archie)



# extract gene sets -----
module_archie <- group_by(signal_archie, Trait, ARCHIE_Component) %>%
  select(Genes) %>%
  mutate(module = cur_group_id()) %>%
  ungroup() %>%
  filter(complete.cases(.))
Nmodule <- n_distinct(module_archie$module)


# extract eQTLGen (trans) summary stat -----
eqtlgen_archie <- eqtlgen[eqtlgen$GeneSymbol %in% module_archie$Genes, ]
eqtlgen_archie$SNPmeta <- paste(eqtlgen_archie$SNPChr, eqtlgen_archie$SNPPos, sep = ":")


## extract common SNPs
res_z <- list()
for(module in 1:Nmodule){
  gene_module <- filter(module_archie, module == !!module) %>% pull(Genes)
  eqtlgen_module <- eqtlgen_archie[eqtlgen_archie$GeneSymbol %in% gene_module, ]
  
  ab_z <- data.table::dcast(
    eqtlgen_module, SNP+SNPmeta~GeneSymbol, fun.aggregate = max,
    value.var = 'Zscore',
    drop = TRUE, fill = NA
  )
  ab_z_na <- is.na(ab_z)
  rem_ind <- rowSums(ab_z_na[, -c(1:2)])
  ab_z_remain <- ab_z[rem_ind == 0, ]
  
  res_z[[str_glue('M{module}')]] <- ab_z_remain
  
  cat("module:", module, "\n")
}



# Indep null SNPs (and #Uniq null SNPs) under various z-score threshold for being null -----
res_nullSNP <- NULL
for(module in 1:Nmodule){
  z_module <- res_z[[str_glue('M{module}')]]
  NSNP <- nrow(z_module); module_size <- ncol(z_module)-2
  
  for(thre_p_z in c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)){
    thre_z2 <- qchisq(1-thre_p_z, 1)
    altSNP_ind <- apply((z_module[, -c(1:2)])^2 > thre_z2, 1, sum) > 0
    
    file_uniq <- paste0("null_SNP/unique.SNP.module", module, ".", thre_p_z, ".txt")
    file_indep <- paste0("null_SNP/indep.SNP.module", module, ".", thre_p_z, ".txt")
    data.table::fwrite(
      z_module[!altSNP_ind, "SNPmeta"], file_uniq,
      quote = FALSE, sep = "\t", col.names = FALSE
    )
    cmd <- str_glue('bash {script_indep_SNP} {file_uniq} {file_indep}')
    system(cmd)
    
    num_nullSNP_indep <- as.numeric(system(str_glue("cat {file_indep} | wc -l"), intern = TRUE))
    res_nullSNP <- rbind(
      res_nullSNP, 
      data.table::data.table(
        "module" = module,
        "module_size" = module_size,
        "thre_z" = thre_p_z,
        "num_nullSNP_indep" = num_nullSNP_indep,
        "num_nullSNP_uniq" = NSNP - sum(altSNP_ind)
      )
    )
  }
}


# print out key message or write out -----
saveRDS(module_archie, file_module_archie)
saveRDS(res_z, file_z)
saveRDS(res_nullSNP, file_num_nullsnp)

