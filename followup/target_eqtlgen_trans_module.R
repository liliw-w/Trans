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




# 1. targets -----
## drug targets
file_gene_repurpose <- '/scratch/midway3/liliw1/chromatin/qa/repurposing_drugs_20200324.txt'

gene_repurpose <- data.table::fread(file_gene_repurpose, skip = 9)

gene_repurpose_expand <- separate_rows(gene_repurpose, target, sep = "[|]")


count(
  gene_repurpose, 
  clinical_phase, 
  if_given_disease_area = disease_area != "", 
  if_corresponding_trait = str_detect(indication, !!indication)
)


filter(
  gene_repurpose_expand, 
  target != "" &
    clinical_phase == "Launched"
) %>%
  distinct(target)

filter(
  gene_repurpose_expand, 
  target != "" &
    clinical_phase == "Launched" &
    disease_area != ""
) %>%
  distinct(target)

filter(
  gene_repurpose_expand, 
  target != "" &
    clinical_phase == "Launched" &
    disease_area != "" &
    str_detect(indication, !!indication)
) %>%
  distinct(target)



dis_cis <- 1e+6/2


# 2. cis snps -----
file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'
file_eqtlgen_snp <- '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'

gene_meta <- data.table::fread(file_gene_meta) %>%
  filter(str_detect(Chromosome, "chr\\d+"))
eqtlgen_snp <- data.table::fread(file_eqtlgen_snp)


target_gene_trait <- filter(
  gene_repurpose_expand, 
  target != "" &
    clinical_phase == "Launched" &
    disease_area != "" &
    str_detect(indication, !!indication)
) %>%
  distinct(target)

sum(target_gene_trait$target %in% gene_meta$GeneSymbol)


target_gene_trait <- inner_join(
  target_gene_trait,
  select(gene_meta, GeneSymbol, Chromosome, Start, End, Class),
  by = c('target' = 'GeneSymbol')
) %>%
  separate(Chromosome, into = c(NA, 'geneChr'), sep = 'chr', convert = TRUE)

target_gene_trait_eqtlgen <- apply(target_gene_trait, 1, function(x){
  tmp_chr = x['geneChr'] |> as.numeric()
  tmp_s = x['Start'] |> as.numeric()
  tmp_e = x['End'] |> as.numeric()
  
  
  filter(eqtlgen_snp, SNPChr == tmp_chr) %>%
    mutate(
      s = abs(SNPPos - tmp_s), 
      e = abs(SNPPos - tmp_e)
    ) %>%
    arrange(s) %>%
    slice(1) %>%
    select(SNPChr, SNPPos, s, e)
}) %>%
  bind_rows() %>%
  cbind(target_gene_trait, .)


sum(target_gene_trait_eqtlgen$s < dis_cis | target_gene_trait_eqtlgen$e < dis_cis)







# 3. coexp modules -----

## if trans snp and trans snp of immune related module -----
file_sig_anno <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/qtl_table.txt'
file_module_meta_enrich <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/module_enrich/module_enrich_top2_immune_related.txt'

sig_anno <- data.table::fread(file_sig_anno)
module_meta_enrich <- data.table::fread(file_module_meta_enrich)


sig_anno <- left_join(
  sig_anno, 
  select(module_meta_enrich, gene_module, if_immune_blood_anno),
  by = c('gene_module')
)

sig_anno <- separate(
  sig_anno, 
  SNP, 
  into = c('SNPChr', 'SNPPos'), 
  sep = ':', 
  remove = FALSE, convert = TRUE
)


target_gene_trait_trans <- apply(target_gene_trait, 1, function(x){
  tmp_chr = x['geneChr'] |> as.numeric()
  tmp_s = x['Start'] |> as.numeric()
  tmp_e = x['End'] |> as.numeric()
  
  tmp_res = select(
    sig_anno,
    gene_module, SNPChr, SNPPos, if_immune_blood_anno
  ) %>%
    filter(SNPChr == tmp_chr) %>%
    mutate(
      s = abs(SNPPos - tmp_s), 
      e = abs(SNPPos - tmp_e)
    ) %>%
    arrange(s) %>%
    filter(
      s < dis_cis | e < dis_cis
    )
  
  tibble(
    "n_trans_module" = pull(tmp_res, gene_module) %>% length(),
    "trans_module" = pull(tmp_res, gene_module) %>% paste(collapse = ";"),
    "n_trans_module_immune" = filter(tmp_res, if_immune_blood_anno) %>% pull(gene_module) %>% length(),
    "trans_module_immune" = filter(tmp_res, if_immune_blood_anno) %>% pull(gene_module) %>% paste(collapse = ";")
  )
}) %>%
  bind_rows() %>%
  cbind(target_gene_trait, .) %>%
  arrange(desc(n_trans_module))


sum(target_gene_trait_trans$n_trans_module > 0)
sum(target_gene_trait_trans$n_trans_module_immune > 0)



data.table::fwrite(
  target_gene_trait_trans,
  file = "target_eqtlgen_trans_module_coexp.txt",
  sep = "\t", quote = FALSE
)





# 3. hallmark modules -----

## if trans snp and trans snp of immune related module -----
file_sig_anno <- '/project2/xuanyao/llw/eQTLGen/MSigDB_est_Sigma/postanalysis/qtl_table.txt'
file_module_meta_enrich <- '/project2/xuanyao/llw/MODULES/MSigDB/msig_50_table_immune_related.txt'


sig_anno <- data.table::fread(file_sig_anno)
module_meta_enrich <- data.table::fread(file_module_meta_enrich)



sig_anno <- left_join(
  sig_anno, 
  select(module_meta_enrich, TERMS, if_immune_blood_anno = if_immune_blood_anno),
  by = c('module_msig' = 'TERMS')
)


sig_anno <- separate(
  sig_anno, 
  SNP, 
  into = c('SNPChr', 'SNPPos'), 
  sep = ':', 
  remove = FALSE, convert = TRUE
)


target_gene_trait_trans2 <- apply(target_gene_trait, 1, function(x){
  tmp_chr = x['geneChr'] |> as.numeric()
  tmp_s = x['Start'] |> as.numeric()
  tmp_e = x['End'] |> as.numeric()
  
  tmp_res = select(
    sig_anno,
    module_msig, SNPChr, SNPPos, if_immune_blood_anno
  ) %>%
    filter(SNPChr == tmp_chr) %>%
    mutate(
      s = abs(SNPPos - tmp_s), 
      e = abs(SNPPos - tmp_e)
    ) %>%
    arrange(s) %>%
    filter(
      s < dis_cis | e < dis_cis
    )
  
  tibble(
    "n_trans_module" = pull(tmp_res, module_msig) %>% length(),
    "trans_module" = pull(tmp_res, module_msig) %>% paste(collapse = ";"),
    "n_trans_module_immune" = filter(tmp_res, if_immune_blood_anno) %>% pull(module_msig) %>% length(),
    "trans_module_immune" = filter(tmp_res, if_immune_blood_anno) %>% pull(module_msig) %>% paste(collapse = ";")
  )
}) %>%
  bind_rows() %>%
  cbind(target_gene_trait, .) %>%
  arrange(desc(n_trans_module))


sum(target_gene_trait_trans2$n_trans_module > 0)
sum(target_gene_trait_trans2$n_trans_module_immune > 0)



data.table::fwrite(
  target_gene_trait_trans2,
  file = "target_eqtlgen_trans_module_hallmark.txt",
  sep = "\t", quote = FALSE
)


# # report -----
# str_glue(
#   "
#   There are {length(unique(gene_repurpose_expand$target))} targets in total.
#   1338 are launched.
#   1153 are with a disease,
#   55 are for the corresponding trait.
#   52 are near a eQTGen snp.
#   
#   30 are associated with at least one coexp-module.
#   21 are associated with at least one immune related coexp-module.
#   
#   20 are associated with at least one hallmark-module.
#   16 are associated with at least one immune related hallmark-module.
#   "
# )

