##############################################
########### Check if drug targets are more likely to be trans signals of gene modules ###########
########### or if trans signals of gene modules are more likely to be drug targets ###########
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

# output
file_res <- str_glue('{trait}_{str_replace(disease_area, "/", "_")}_{indication}_target_enrich.txt')



## 1. select snp -----
gwas <- data.table::fread(
  file = file_gwas,
  select = gwas_col,
  col.names = c('CHR', 'POS', 'gwas_pval')
)
gwas$SNP <- paste(gwas$CHR, gwas$POS, sep = ":")


# in eQTLGen
file_eqtlgen_snp <- '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'
eqtlgen_snp <- data.table::fread(file_eqtlgen_snp)

gwas_eqtlgen <- filter(
  gwas,
  SNP %in% !!eqtlgen_snp$meta
)


# trait significant
gwas_eqtlgen$if_trait_sig <- gwas_eqtlgen$gwas_pval < 0.05/(1e+6)
gwas_eqtlgen_sig <- filter(gwas_eqtlgen, if_trait_sig)




## 2. add nearest gene -----
file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
dis_cis <- 1e+6

gene.meta <- data.table::fread(file_gene_meta, header = TRUE)

gene.meta <- gene.meta %>%
  filter(Class %in% c("protein_coding", "lincRNA") & Chromosome %in% paste0("chr",1:22)) %>%
  separate("Chromosome", c(NA, "chr"), sep = "chr", remove = FALSE)

cis.gene.meta <- sapply(1:nrow(gwas_eqtlgen_sig), function(x){
  tmp_qtl = gwas_eqtlgen_sig[x, ]
  tmp_cis.gene.meta = gene.meta %>%
    filter(chr %in% tmp_qtl$CHR) %>%
    mutate("dis" = abs(tmp_qtl$POS - Start)) %>%
    filter(dis < dis_cis/2) %>%
    arrange(dis)
  c(
    tmp_cis.gene.meta$GeneSymbol[1],
    tmp_cis.gene.meta$dis[1],
    paste(tmp_cis.gene.meta$GeneSymbol, collapse = ";"),
    paste(tmp_cis.gene.meta$dis, collapse = ";")
  )
})
cis.gene.meta <- data.table::as.data.table(t(cis.gene.meta))
colnames(cis.gene.meta) <- c("nearest_gene", "nearest_dis", "near_genes", "near_dis")

## add cis gene info to signals
gwas_eqtlgen_sig <- bind_cols(gwas_eqtlgen_sig, cis.gene.meta)



## 3. if nearest gene is target -----
file_gene_repurpose <- '/scratch/midway3/liliw1/chromatin/qa/repurposing_drugs_20200324.txt'

gene_repurpose <- data.table::fread(file_gene_repurpose, skip = 9)

target_gene_all <- filter(gene_repurpose, disease_area != "") %>%
  filter(target != "") %>%
  pull(target) %>%
  str_split(., pattern = "\\|") %>%
  unlist() %>%
  unique()

target_gene_trait <- filter(gene_repurpose, str_detect(indication, !!indication) ) %>%
  filter(target != "") %>%
  pull(target) %>%
  str_split(., pattern = "\\|") %>%
  unlist() %>%
  unique()


gwas_eqtlgen_sig$if_target_all <- gwas_eqtlgen_sig$nearest_gene %in% target_gene_all
gwas_eqtlgen_sig$if_target_trait <- gwas_eqtlgen_sig$nearest_gene %in% target_gene_trait



## 4. if trans snp (2161) and trans snp of immune related module -----
file_sig_anno <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/qtl_table.txt'
file_module_meta_enrich <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich_top2_immune_related.txt'

sig_anno <- data.table::fread(file_sig_anno)
module_meta_enrich <- data.table::fread(file_module_meta_enrich)


sig_anno <- left_join(
  sig_anno, 
  select(module_meta_enrich, gene_module, if_immune_blood_anno),
  by = c('gene_module')
)

# add if trans snp
trans_unique_snp <- unique(sig_anno$SNP)
gwas_eqtlgen_sig$if_trans_sig <- gwas_eqtlgen_sig$SNP %in% trans_unique_snp$snp


# add if trans snp of immune related module
# to limit the analyses to trans-eQTLs of immune related modules only
gwas_eqtlgen_sig <- left_join(
  gwas_eqtlgen_sig, 
  
  filter(sig_anno, if_immune_blood_anno) %>%
    group_by(SNP) %>%
    summarise(
      if_immune_blood_anno = sum(if_immune_blood_anno),
      gene_module_immune = paste(gene_module, collapse = ";")
    ) %>%
    ungroup(),
  
  by = c('SNP')
) %>%
  mutate(
    if_immune_blood_anno = replace(if_immune_blood_anno, is.na(if_immune_blood_anno), FALSE),
    gene_module_immune = replace(gene_module_immune, is.na(gene_module_immune), "")
  )



## 5. enrichment -----
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
  mutate(gwas_eqtlgen_sig, "trait" = !!trait),
  file = file_res,
  sep = "\t", quote = FALSE
)

