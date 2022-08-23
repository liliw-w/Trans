rm(list = ls())
library(data.table)
library(tidyverse)


#file_pheno_manifest = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"
#dir_gwas_data = "/project2/xuanyao/llw/GWAS/UKB_nealelab"
#dir_coloc_data = "/project2/xuanyao/llw/coloc"
#dir_coloc = "/scratch/midway2/liliw1/coloc"

#pheno_manifest = fread(file_pheno_manifest)
#pheno_manifest$filename

########## files and parameters ##########
pop = "EUR"
gwasPhenocode = ""
gwas_trait_type = ""
gwas_trait = "height_Meta-analysis_Wood_et_al+UKBiobank_2018"

#gwasPhenocode_seq = pheno_manifest$phenocode_uniq
p_included_thre = 1e-5

file_qtlColocReg = "/project2/xuanyao/llw/coloc/qtlColocReg.txt.gz"
#file_snp_meta = file.path(dir_gwas_data, "full_variant_qc_metrics.txt.bgz")
file_gwas_bgz = "/project2/xuanyao/llw/GWAS/height_Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz"


#dir_gwas = file.path(dir_coloc,
#                     paste("ukbb",
#                           gwas_trait_type,
#                           gwasPhenocode,
#                           gwas_trait,
#                           sep = "_")
#)
#dir.create(file.path(dir_gwas, "data"), recursive = TRUE, showWarnings = FALSE)


## output files
#setwd(dir_gwas)
file_qtlColocReg_gwas = "qtlColocReg_gwas.txt.gz"
file_gwasColocReg = "gwasColocReg.txt.gz"
file_gwasRegTruncPthre = "gwasRegTruncPthre.txt"




########## read files ##########
#gwasCol = c("chr", "pos", paste(c("af", "beta", "se", "pval", "low_confidence"), pop, sep = "_") )
gwas = fread(cmd = paste("gunzip -c", file_gwas_bgz))

#snpMetaCol = c("chrom", "pos", "rsid", "high_quality")
#snpMeta = fread(cmd = paste("gunzip -c", file_snp_meta), select = snpMetaCol)

qtlColocReg = fread(file_qtlColocReg)


########## SNP filtering for gwas ##########
#if(identical(gwas$pos, snpMeta$pos)){
#  gwas = cbind(gwas, snpMeta[, c("rsid", "high_quality")])
#}else stop("GWAS variants and SNP meta are not aligned!")


#gwas$if_good = !is.na(gwas[[paste0("pval_", pop)]]) & gwas$high_quality & !gwas[[paste0("low_confidence_", pop)]]


########## gwas good SNPs & if overlapped with qtl SNPs & remove duplicated SNPs ##########
gwas$if_qtl = paste(gwas$CHR, gwas$POS, sep = ":") %in% qtlColocReg$SNP_ID
gwas = gwas %>% filter(if_qtl &
                         !duplicated(paste(CHR, POS, sep = ":")) )
gwas$SNP_ID = paste(gwas$CHR, gwas$POS, sep = ":")

if(nrow(gwas) == 0) stop("No individuals in the specified population are phenotyped for this GWAS trait!")


########## qtl regions with overlapped SNPs ##########
qtlColocReg_gwas = qtlColocReg[qtlColocReg$SNP_ID %in% gwas$SNP_ID, ]


########## coloc regions based on qtl ##########
#gwas$chr = as.integer(gwas$chr)

gwasColocReg = qtlColocReg_gwas %>% select(c(Signal, Module, SNP_ID, Chr, Pos, Region)) %>% left_join(y = gwas, by = c("SNP_ID", "Chr"="CHR", "Pos"="POS"))


# range of pvalues in a region for GWAS and qtl & remove the regions whose min GWAS p < 1e-5
gwasRegP = gwasColocReg %>% group_by(Region) %>% summarise(s = min(P), l = max(P))
qtlRegP = qtlColocReg_gwas %>% group_by(Region) %>% summarise(s = min(Pval), l = max(Pval))

gwasRegTruncPthre = gwasRegP$Region[gwasRegP$s < p_included_thre]


########## save results ##########
fwrite(qtlColocReg_gwas, file_qtlColocReg_gwas, quote = FALSE, sep = "\t")
fwrite(gwasColocReg, file_gwasColocReg, quote = FALSE, sep = "\t")
fwrite(data.table(gwasRegTruncPthre), file_gwasRegTruncPthre, quote = FALSE, sep = "\t", col.names = FALSE)


### Notes

# for trait 30080 Platelet count, there are 28,987,535 rows.
# For snp meta, there are 28,987,535 rows.

### only considering high-quality variants, `high_quality`
# out of 28987534 snps, 26691267 (~92%) are high_quality.

### retained variants with INFO scores > 0.8 (29,865,259 variants on the autosomes and X chromosome)
# For all snps in snp_meta file, the minimum info score is 0.8. So these snps are already filtered.

### low-confidence statistics removed, `low_confidence_{pop}`
