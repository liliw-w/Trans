##############################################
########### Run trans-pco on archie gene sets ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)


# I/O & paras -----
thre_p_z <- 1e-4
dir_pco <- '/scratch/midway3/liliw1/chromatin/script/pco/'

file_module_archie <- 'null_SNP/module_archie.rds'
file_eQTLGen_z <- "null_SNP/sumStat_archie.rds"
file_null_SNP <- "null_SNP/num_nullSNP.rds"
file_gene_meta <- '/project2/xuanyao/data/eQTLGen/gene_meta_trans.txt.gz'

file_indep_null_snp <- list.files(
  'null_SNP',
  str_glue('indep\\.SNP\\.module\\d+\\.{thre_p_z}\\.txt'),
  full.names = TRUE
)

## output -----
file_p_all <- str_glue('null_SNP/p_all_archie.txt.gz')



# read files -----
module_archie <- readRDS(file_module_archie)
eQTLGen_z <- readRDS(file_eQTLGen_z)
res_nullSNP <- readRDS(file_null_SNP)
gene_meta <- data.table::fread(file_gene_meta, header = TRUE)
indep_null_snp <- lapply(file_indep_null_snp, data.table::fread, header = FALSE, col.names = "snp") %>%
  setNames(
    str_extract(file_indep_null_snp, 'module\\d+') %>% str_extract('\\d+') %>% as.numeric(.)
  ) %>%
  bind_rows(.id = "module")


source(paste0(dir_pco, "ModifiedPCOMerged_acat.R"))
source(paste0(dir_pco, "liu.R"))
source(paste0(dir_pco, "liumod.R"))
source(paste0(dir_pco, "davies.R"))
dyn.load(paste0(dir_pco, "qfc.so"))



# check how many snps will be used for this module under given null p threshold -----
res_nullSNP %>% filter(thre_z == thre_p_z)


# trans-pco on archie gene sets -----
Nmodule <- n_distinct(module_archie$module)
p_all <- list()
for(module in 1:Nmodule){
  eQTLGen_z_module <- eQTLGen_z[[str_glue('M{module}')]]
  indep_null_snp_module <- filter(indep_null_snp, module == !!module)
  gene_meta_module <- filter(gene_meta, GeneSymbol %in% colnames(eQTLGen_z_module))
  
  
  ## Estimate sigma using extracted independent null z for each module -----
  indep_null_z <- eQTLGen_z_module[eQTLGen_z_module$SNPmeta %in% indep_null_snp_module$snp, -c(1:2)]
  Sigma_null_z <- cor(indep_null_z)
  
  
  ## extracted only trans variants either more 5Mb away from all genes or on a different chromosome -----
  eQTLGen_z_module$if_trans <- select(eQTLGen_z_module, SNPmeta) %>%
    separate(
      SNPmeta, into = c('snp_chr', 'snp_pos'), sep = ":", remove = FALSE
    ) %>%
    apply(., 1, function(x){
      filter(
        gene_meta_module, 
        GeneChr == as.numeric(x['snp_chr']) & 
          abs(GenePos - as.numeric(x['snp_pos'])) < 5e+6
      ) %>%
        nrow(.) == 0
    })
  
  z_mat <- filter(eQTLGen_z_module, if_trans) %>%
    select(!c(SNP, SNPmeta, if_trans)) %>%
    as.matrix()
  rownames(z_mat) <- eQTLGen_z_module$SNP
  
  
  ## pco -----
  p_all[[str_glue('M{module}')]] <- ModifiedPCOMerged_acat(Z.mat = z_mat, Sigma = Sigma_null_z) %>%
    enframe("snp", "p")
  
  cat(str_glue('M{module} done out of {Nmodule} modules. \n\n'))
}


## add module meta -----
p_all <- bind_rows(p_all, .id = "module") %>%
  left_join(
    distinct(module_archie, Trait, ARCHIE_Component, module = str_glue('M{module}')),
    by = c('module')
  )



# print out key message or write out -----
data.table::fwrite(
  p_all,
  file_p_all,
  quote = FALSE, sep = '\t'
)


# # Vis #independent null SNPs for estimating $\Sigma_{nullZ}$ v.s. gene set size -----
# left_join(
#   count(module_archie, Trait, ARCHIE_Component, module, name = "archie_gene_set"),
#   count(indep_null_snp, module = as.numeric(module), name = "indep_null_snps"),
#   by = c('module')
# ) %>%
#   pivot_longer(
#     cols = c(archie_gene_set, indep_null_snps),
#     names_to = "Type", values_to = "num_snp"
#   ) %>%
#   ggplot(aes(x = str_glue('{Trait}_{ARCHIE_Component}'), y = num_snp, fill = Type)) +
#   geom_col(position = "dodge") +
#   geom_text(aes(label = num_snp), position = position_dodge(width = .9)) +
#   scale_fill_brewer(palette = "Paired") +
#   labs(x = "Archie signal set", y = "#genes/SNPs") +
#   theme(
#     axis.text.x = element_text(angle = 90)
#   )
# 
# ggsave(
#   filename = 'null_SNP/plt_gene_indepnull_size.pdf',
#   height = 3, width = 5
# )
# 

