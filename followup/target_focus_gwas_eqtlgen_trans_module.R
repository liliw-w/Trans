rm(list = ls())
library(tidyverse)

dis_cis <- 1e+6/2

k <- 5

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



# 0. eqtlgen 10k snps -----
file_eqtlgen_snp <- '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'

eqtlgen_snp <- data.table::fread(file_eqtlgen_snp)


# 1. if target genes -----
## drug targets
file_gene_repurpose <- '/scratch/midway3/liliw1/chromatin/qa/repurposing_drugs_20200324.txt'
file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'

gene_repurpose <- data.table::fread(file_gene_repurpose, skip = 9)
gene_meta <- data.table::fread(file_gene_meta) %>%
  filter(str_detect(Chromosome, "chr\\d+")) %>%
  filter(Class %in% c('protein_coding', 'lincRNA'))


gene_repurpose_expand <- separate_rows(gene_repurpose, target, sep = "[|]")
gene_meta <- select(gene_meta, GeneSymbol, Chromosome, Start, End, Class) %>%
  separate(Chromosome, into = c(NA, 'geneChr'), sep = 'chr', convert = TRUE)


target_gene_trait <- filter(
  gene_repurpose_expand, 
  target != "" &
    clinical_phase == "Launched" &
    disease_area != "" &
    str_detect(indication, !!indication)
) %>%
  distinct(target)

sum(target_gene_trait$target %in% gene_meta$GeneSymbol)
sum(gene_meta$GeneSymbol %in% target_gene_trait$target)


gene_meta$if_target <- gene_meta$GeneSymbol %in% target_gene_trait$target



target_gene_all <- filter(
  gene_repurpose_expand, 
  target != "" &
    clinical_phase == "Launched" &
    disease_area != ""
) %>%
  distinct(target)

gene_meta$if_target_any <- gene_meta$GeneSymbol %in% target_gene_all$target





# 2. cis genes -----
eqtlgen_snp_res <- apply(eqtlgen_snp, 1, function(x){
  tmp_chr = x['SNPChr'] |> as.numeric()
  tmp_pos = x['SNPPos'] |> as.numeric()
  
  tmp_res = filter(gene_meta, geneChr == tmp_chr) %>%
    mutate(
      s = abs(tmp_pos - Start), 
      e = abs(tmp_pos - End),
      minse = min(c(s, e)),
    ) %>%
    arrange(s) %>%
    filter(
      s < dis_cis | e < dis_cis
    )
  
  
  tibble(
    "n_cis_gene" = pull(tmp_res, GeneSymbol) %>% length(),
    "cis_gene" = pull(tmp_res, GeneSymbol) %>% paste(collapse = ";"),
    "n_target" = sum(tmp_res$if_target),
    "if_target" = pull(tmp_res, if_target) %>% paste(collapse = ";"),
    "n_target_any" = sum(tmp_res$if_target_any),
    "if_target_any" = pull(tmp_res, if_target_any) %>% paste(collapse = ";")
  )
}) %>%
  bind_rows() %>%
  cbind(eqtlgen_snp, .)



# 3. if trait associated snps -----
gwas <- data.table::fread(
  file = file_gwas,
  select = gwas_col,
  col.names = c('CHR', 'POS', 'gwas_pval')
)

# trait significant
gwas$SNP <- paste(gwas$CHR, gwas$POS, sep = ":")
gwas_sig <- filter(gwas, gwas_pval < 0.05/(1e+6))

eqtlgen_snp_res$if_trait_sig <- eqtlgen_snp_res$meta %in% gwas_sig$SNP




# 4. if trans snps associated with immune modules -----

## 4.1. co-exp modules -----
file_sig_anno <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/qtl_table.txt'
file_module_meta_enrich <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich_top2_immune_related.txt'

sig_anno <- data.table::fread(file_sig_anno)
module_meta_enrich <- data.table::fread(file_module_meta_enrich)

sig_anno <- left_join(
  sig_anno, 
  select(module_meta_enrich, gene_module, if_immune_blood_anno),
  by = c('gene_module')
)


eqtlgen_snp_res$if_trans_immune_coexp <- eqtlgen_snp_res$meta %in% (filter(sig_anno, if_immune_blood_anno) %>% distinct(SNP) %>% pull(SNP))
eqtlgen_snp_res$if_trans_sig_coexp <- eqtlgen_snp_res$meta %in% (distinct(sig_anno, SNP) %>% pull(SNP))

eqtlgen_snp_res <- left_join(
  eqtlgen_snp_res,
  
  select(
    sig_anno,
    gene_module, SNP
  ) %>%
    group_by(SNP) %>%
    summarise(module_sig_coexp = paste(gene_module, collapse = "_")) %>%
    ungroup(),
  
  by = c('meta' = 'SNP')
) %>%
  left_join(
    filter(sig_anno, if_immune_blood_anno) %>%
      select(
        gene_module, SNP
      ) %>%
      group_by(SNP) %>%
      summarise(module_immune_coexp = paste(gene_module, collapse = "_")) %>%
      ungroup(),
    
    by = c('meta' = 'SNP')
  )



## 4.2. hallmark modules -----
file_sig_anno <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/qtl_table.txt'
file_module_meta_enrich <- '/project2/xuanyao/llw/MODULES/MSigDB/msig_50_table_immune_related.txt'

sig_anno <- data.table::fread(file_sig_anno)
module_meta_enrich <- data.table::fread(file_module_meta_enrich)

sig_anno <- left_join(
  sig_anno, 
  select(module_meta_enrich, TERMS, if_immune_blood_anno = if_immune_blood_anno),
  by = c('module_msig' = 'TERMS')
)


eqtlgen_snp_res$if_trans_immune_hallmark <- eqtlgen_snp_res$meta %in% (filter(sig_anno, if_immune_blood_anno) %>% distinct(SNP) %>% pull(SNP))
eqtlgen_snp_res$if_trans_sig_hallmark <- eqtlgen_snp_res$meta %in% (distinct(sig_anno, SNP) %>% pull(SNP))


eqtlgen_snp_res <- left_join(
  eqtlgen_snp_res,
  
  select(
    sig_anno,
    module_msig, SNP
  ) %>%
    group_by(SNP) %>%
    summarise(module_sig_hallmark = paste(module_msig, collapse = "_")) %>%
    ungroup(),
  
  by = c('meta' = 'SNP')
) %>%
  left_join(
    filter(sig_anno, if_immune_blood_anno) %>%
      select(
        module_msig, SNP
      ) %>%
      group_by(SNP) %>%
      summarise(module_immune_hallmark = paste(module_msig, collapse = "_")) %>%
      ungroup(),
    
    by = c('meta' = 'SNP')
  )


# data.table::fwrite(
#   eqtlgen_snp_res,
#   file = "eqtlgen_snp_target_gwas_trans_module.txt",
#   sep = "\t", quote = FALSE
# )




# 5. enrich (gene based) -----
count(eqtlgen_snp_res, n_target)
count(eqtlgen_snp_res, if_trait_sig)
count(eqtlgen_snp_res, if_trans_immune_coexp)
count(eqtlgen_snp_res, if_trans_immune_hallmark)


target_focus_res <- separate_rows(eqtlgen_snp_res, cis_gene, if_target, sep = "[;]", convert = TRUE) %>%
  filter(cis_gene != "") %>%
  group_by(cis_gene, if_target) %>%
  summarise(
    SNP = paste(meta, collapse = ";"),
    module_immune_coexp = paste(unique(module_immune_coexp), collapse = ";"),
    module_immune_hallmark = paste(unique(module_immune_hallmark), collapse = ";"),
    if_trans_sig_coexp = any(if_trans_sig_coexp),
    if_trans_sig_hallmark = any(if_trans_sig_hallmark),
    if_trait_sig = any(if_trait_sig),
    if_trans_immune_coexp = any(if_trans_immune_coexp),
    if_trans_immune_hallmark = any(if_trans_immune_hallmark)
  ) %>%
  ungroup()
target_focus_res$if_trans_immune_any <- target_focus_res$if_trans_immune_coexp | target_focus_res$if_trans_immune_hallmark



## enrich coexp
enrich_mat_coexp <- matrix(
  c(
    sum(target_focus_res$if_trans_immune_coexp & target_focus_res$if_trait_sig & target_focus_res$if_target),
    sum(!target_focus_res$if_trans_immune_coexp & target_focus_res$if_trait_sig & target_focus_res$if_target),
    sum(target_focus_res$if_trans_immune_coexp & (!target_focus_res$if_trait_sig & !target_focus_res$if_target)),
    sum(!target_focus_res$if_trans_immune_coexp & (!target_focus_res$if_trait_sig & !target_focus_res$if_target))
  ),
  ncol = 2
)

enrich_mat_coexp
fisher.test(enrich_mat_coexp, alternative = "greater")

filter(
  target_focus_res,
  if_trait_sig & if_target
) %>%
  View()


## enrich hallmark
enrich_mat_hallmark <- matrix(
  c(
    sum(target_focus_res$if_trans_immune_hallmark & target_focus_res$if_trait_sig & target_focus_res$if_target),
    sum(!target_focus_res$if_trans_immune_hallmark & target_focus_res$if_trait_sig & target_focus_res$if_target),
    sum(target_focus_res$if_trans_immune_hallmark & (!target_focus_res$if_trait_sig & !target_focus_res$if_target)),
    sum(!target_focus_res$if_trans_immune_hallmark & (!target_focus_res$if_trait_sig & !target_focus_res$if_target))
  ),
  ncol = 2
)

enrich_mat_hallmark
fisher.test(enrich_mat_hallmark, alternative = "greater")


## enrich any
enrich_mat_any <- matrix(
  c(
    sum(target_focus_res$if_trans_immune_any & target_focus_res$if_trait_sig & target_focus_res$if_target),
    sum(!target_focus_res$if_trans_immune_any & target_focus_res$if_trait_sig & target_focus_res$if_target),
    sum(target_focus_res$if_trans_immune_any & (!target_focus_res$if_trait_sig & !target_focus_res$if_target)),
    sum(!target_focus_res$if_trans_immune_any & (!target_focus_res$if_trait_sig & !target_focus_res$if_target))
  ),
  ncol = 2
)

enrich_mat_any
fisher.test(enrich_mat_any, alternative = "greater")



# 6. enrich (snp based) -----
snp_focus_res <- filter(eqtlgen_snp_res, cis_gene != "")

snp_focus_res$if_trans_immune_any <- snp_focus_res$if_trans_immune_coexp | snp_focus_res$if_trans_immune_hallmark


## enrich coexp
enrich_mat_coexp <- matrix(
  c(
    sum(snp_focus_res$if_trans_immune_coexp & snp_focus_res$if_trait_sig & snp_focus_res$n_target),
    sum(!snp_focus_res$if_trans_immune_coexp & snp_focus_res$if_trait_sig & snp_focus_res$n_target),
    sum(snp_focus_res$if_trans_immune_coexp & (!snp_focus_res$if_trait_sig & !snp_focus_res$n_target)),
    sum(!snp_focus_res$if_trans_immune_coexp & (!snp_focus_res$if_trait_sig & !snp_focus_res$n_target))
  ),
  ncol = 2
)

enrich_mat_coexp
fisher.test(enrich_mat_coexp, alternative = "greater")



## enrich hallmark
enrich_mat_hallmark <- matrix(
  c(
    sum(snp_focus_res$if_trans_immune_hallmark & snp_focus_res$if_trait_sig & snp_focus_res$n_target),
    sum(!snp_focus_res$if_trans_immune_hallmark & snp_focus_res$if_trait_sig & snp_focus_res$n_target),
    sum(snp_focus_res$if_trans_immune_hallmark & (!snp_focus_res$if_trait_sig & !snp_focus_res$n_target)),
    sum(!snp_focus_res$if_trans_immune_hallmark & (!snp_focus_res$if_trait_sig & !snp_focus_res$n_target))
  ),
  ncol = 2
)

enrich_mat_hallmark
fisher.test(enrich_mat_hallmark, alternative = "greater")


## enrich any
enrich_mat_any <- matrix(
  c(
    sum(snp_focus_res$if_trans_immune_any & snp_focus_res$if_trait_sig & snp_focus_res$n_target),
    sum(!snp_focus_res$if_trans_immune_any & snp_focus_res$if_trait_sig & snp_focus_res$n_target),
    sum(snp_focus_res$if_trans_immune_any & (!snp_focus_res$if_trait_sig & !snp_focus_res$n_target)),
    sum(!snp_focus_res$if_trans_immune_any & (!snp_focus_res$if_trait_sig & !snp_focus_res$n_target))
  ),
  ncol = 2
)

enrich_mat_any
fisher.test(enrich_mat_any, alternative = "greater")




# 7. trans sig based -----
## 1. enrich (gene based) -----
target_focus_res <- separate_rows(eqtlgen_snp_res, cis_gene, if_target, if_target_any, sep = "[;]", convert = TRUE) %>%
  filter(cis_gene != "") %>%
  group_by(cis_gene, if_target, if_target_any) %>%
  summarise(
    SNP = paste(meta, collapse = ";"),
    module_immune_coexp = paste(unique(module_immune_coexp), collapse = ";"),
    module_immune_hallmark = paste(unique(module_immune_hallmark), collapse = ";"),
    if_trait_sig = any(if_trait_sig),
    if_trans_immune_coexp = any(if_trans_immune_coexp),
    if_trans_immune_hallmark = any(if_trans_immune_hallmark),
    if_trans_sig_coexp = any(if_trans_sig_coexp),
    if_trans_sig_hallmark = any(if_trans_sig_hallmark)
  ) %>%
  ungroup()
target_focus_res$if_trans_immune_any <- target_focus_res$if_trans_immune_coexp | target_focus_res$if_trans_immune_hallmark
target_focus_res$if_trans_sig_any <- target_focus_res$if_trans_sig_coexp | target_focus_res$if_trans_sig_hallmark


sum(target_focus_res$if_target & target_focus_res$if_trait_sig)

filter(
  target_focus_res,
  if_trans_sig_any & if_trans_immune_any & if_target & if_trait_sig
) %>%
  select(cis_gene, module_immune_coexp, module_immune_hallmark, SNP)



## enrich any
enrich_mat_any <- matrix(
  c(
    sum(target_focus_res$if_trans_sig_any & target_focus_res$if_trans_immune_any & target_focus_res$if_target & target_focus_res$if_trait_sig),
    sum(target_focus_res$if_trans_sig_any & !target_focus_res$if_trans_immune_any & target_focus_res$if_target & target_focus_res$if_trait_sig),
    sum(target_focus_res$if_trans_sig_any & target_focus_res$if_trans_immune_any & !target_focus_res$if_target & !target_focus_res$if_trait_sig),
    sum(target_focus_res$if_trans_sig_any & !target_focus_res$if_trans_immune_any & !target_focus_res$if_target & !target_focus_res$if_trait_sig)
  ),
  ncol = 2
)

enrich_mat_any
fisher.test(enrich_mat_any, alternative = "greater")



## 2. enrich (snp based) -----
snp_focus_res <- filter(eqtlgen_snp_res, if_trans_sig_coexp | if_trans_sig_hallmark) %>%
  filter(cis_gene != "")

snp_focus_res$if_trans_immune_any <- snp_focus_res$if_trans_immune_coexp | snp_focus_res$if_trans_immune_hallmark


## enrich any
enrich_mat_any <- matrix(
  c(
    sum(snp_focus_res$if_trans_immune_any & snp_focus_res$n_target & snp_focus_res$n_target_any),
    sum(!snp_focus_res$if_trans_immune_any & snp_focus_res$n_target & snp_focus_res$n_target_any),
    sum(snp_focus_res$if_trans_immune_any & !snp_focus_res$n_target & !snp_focus_res$n_target_any),
    sum(!snp_focus_res$if_trans_immune_any & !snp_focus_res$n_target & !snp_focus_res$n_target_any)
  ),,
  ncol = 2
)

enrich_mat_any
fisher.test(enrich_mat_any, alternative = "greater")



