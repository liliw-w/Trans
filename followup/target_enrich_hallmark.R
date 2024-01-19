##############################################
########### Check if drug targets are more likely to be trans signals of gene modules ###########
########### or if trans signals of gene modules are more likely to be drug targets ###########
########### use hallmark gene sets as modules ###########
##############################################
rm(list = ls())
library(tidyverse)


k <- 1

# trait meta
all_trait_meta <- tibble(
  "trait" = c('allergy', "cd", 'ibd', "UC", "MS"),
  "disease_area" = c('allergy', "gastroenterology", "gastroenterology", "gastroenterology", "neurology/psychiatry"),
  "indication" = c('allergic', 'Crohn', 'inflammatory bowel disease', 'ulcerative colitis', 'multiple sclerosis'),
  "file_gwas" = c(
    '/project2/xuanyao/llw/GWAS/pmid29083406_Allergy_SHARE-without23andMe.LDSCORE-GC.SE-META.v0_add_s.gz',
    '/project2/xuanyao/llw/GWAS/pmid28067908_cd_build37_40266_20161107.txt.gz',
    '/project2/xuanyao/llw/GWAS/pmid28067908_ibd_build37_59957_20161107.txt.gz',
    '/project2/xuanyao/llw/GWAS/pmid26192919_EUR.UC.gwas_info03_filtered.assoc.gz',
    '/project2/xuanyao/llw/GWAS/pmid31604244_MS_15_discovery_metav3.0.meta.gz'
  ),
  'gwas_col' = list(
    c('CHR', 'BP', 'PVALUE'),
    c('Chr', 'Pos', 'P.value'),
    c('Chr', 'Pos', 'P.value'),
    c('CHR', 'BP', 'P'),
    c('CHR', 'BP', 'P')
  )
)

trait <- all_trait_meta$trait[k]
disease_area <- all_trait_meta$disease_area[k]
indication <- all_trait_meta$indication[k]
file_gwas <- all_trait_meta$file_gwas[k]
gwas_col <- all_trait_meta$gwas_col[[k]]

file_original <- str_glue('{trait}_{str_replace(disease_area, "/", "_")}_{indication}_target_enrich.txt')

# output
file_res <- str_glue('{trait}_{str_replace(disease_area, "/", "_")}_{indication}_target_enrich_hallmark.txt')


## 1. old res file -----
gwas_eqtlgen_sig <- data.table::fread(file_original)



## 2. if trans snp and trans snp of immune related module (hallmark) -----
file_sig_anno <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/qtl_table.txt'
file_module_meta_enrich <- '/project2/xuanyao/llw/MODULES/MSigDB/msig_50_table_immune_related.txt'


sig_anno <- data.table::fread(file_sig_anno)
module_meta_enrich <- data.table::fread(file_module_meta_enrich)


sig_anno <- left_join(
  sig_anno, 
  select(module_meta_enrich, TERMS, if_immune_blood_anno = if_immune_blood_anno),
  by = c('module_msig' = 'TERMS')
)

# add if trans snp
trans_unique_snp <- unique(sig_anno$SNP)
gwas_eqtlgen_sig$if_trans_sig <- gwas_eqtlgen_sig$SNP %in% trans_unique_snp


# add if trans snp of immune related module
# to limit the analyses to trans-eQTLs of immune related modules only
gwas_eqtlgen_sig <- left_join(
  select(
    gwas_eqtlgen_sig,
    !c(if_immune_blood_anno, gene_module_immune)
  ), 
  
  filter(sig_anno, if_immune_blood_anno) %>%
    group_by(SNP) %>%
    summarise(
      if_immune_blood_anno = sum(if_immune_blood_anno),
      gene_module_immune = paste(module_msig, collapse = ";")
    ) %>%
    ungroup(),
  
  by = c('SNP')
) %>%
  mutate(
    if_immune_blood_anno = replace(if_immune_blood_anno, is.na(if_immune_blood_anno), FALSE),
    gene_module_immune = replace(gene_module_immune, is.na(gene_module_immune), "")
  )



## 3. enrichment -----
# use all targets
enrich_mat_target_all <- count(gwas_eqtlgen_sig, if_trans_sig, if_target_all) %>%
  arrange(desc(if_trans_sig), desc(if_target_all))
enrich_mat_target_all

if(nrow(enrich_mat_target_all) == 4) fisher.test(matrix(enrich_mat_target_all$n, ncol = 2), alternative = "greater")


# use trait specific targets
enrich_mat_target_trait <- count(gwas_eqtlgen_sig, if_trans_sig, if_target_trait) %>%
  arrange(desc(if_trans_sig), desc(if_target_trait))
enrich_mat_target_trait

if(nrow(enrich_mat_target_trait) == 4) fisher.test(matrix(enrich_mat_target_trait$n, ncol = 2), alternative = "greater")



# save -----
data.table::fwrite(
  gwas_eqtlgen_sig,
  file = file_res,
  sep = "\t", quote = FALSE
)

