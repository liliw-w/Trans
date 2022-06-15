rm(list = ls())
library(data.table)
library(tidyverse)


file_pheno_manifest = "/scratch/midway2/liliw1/coloc/phenotype_manifest_sub.tsv"
dir_gwas_data = "/project2/xuanyao/llw/GWAS/UKB_nealelab"
dir_coloc_data = "/project2/xuanyao/llw/coloc"
dir_coloc = "/scratch/midway2/liliw1/coloc"

pheno_manifest = fread(file_pheno_manifest)
pheno_manifest$filename

########## files and parameters ##########
pop = "EUR"
gwasPhenocode = pheno_manifest$phenocode_uniq[1]
gwas_trait_type = pheno_manifest$trait_type[1]
gwas_trait = pheno_manifest$trait[1]

#gwasPhenocode_seq = pheno_manifest$phenocode_uniq
p_included_thre = 1e-5

file_qtlColocReg = file.path(dir_coloc_data, "qtlColocReg.txt.gz")
file_snp_meta = file.path(dir_gwas_data, "full_variant_qc_metrics.txt.bgz")
file_gwas_bgz = file.path(dir_gwas_data, pheno_manifest$filename)[1]


dir_gwas = file.path(dir_coloc,
                     paste("ukbb",
                           gwas_trait_type,
                           gwasPhenocode,
                           gwas_trait,
                           sep = "_")
                     )
dir.create(file.path(dir_gwas, "data"), recursive = TRUE, showWarnings = FALSE)


## output files
setwd(dir_gwas)
file_qtlColocReg_gwas = "data/qtlColocReg_gwas.txt.gz"
file_gwasColocReg = "data/gwasColocReg.txt.gz"
file_gwasRegTruncPthre = "data/gwasRegTruncPthre.txt"




########## read files ##########
gwasCol = c("chr", "pos", paste(c("af", "beta", "se", "pval", "low_confidence"), pop, sep = "_") )
gwas = fread(cmd = paste("gunzip -c", file_gwas_bgz), select = gwasCol)

snpMetaCol = c("chrom", "pos", "rsid", "high_quality")
snpMeta = fread(cmd = paste("gunzip -c", file_snp_meta), select = snpMetaCol)

qtlColocReg = fread(file_qtlColocReg)


########## SNP filtering for gwas ##########
if(identical(gwas$pos, snpMeta$pos)){
  gwas = cbind(gwas, snpMeta[, c("rsid", "high_quality")])
}else stop("GWAS variants and SNP meta are not aligned!")


gwas$if_good = !is.na(gwas[[paste0("pval_", pop)]]) & gwas$high_quality & !gwas[[paste0("low_confidence_", pop)]]


########## gwas good SNPs & if overlapped with qtl SNPs & remove duplicated SNPs ##########
gwas$if_qtl = paste(gwas$chr, gwas$pos, sep = ":") %in% qtlColocReg$SNP_ID
gwas = gwas %>% filter(if_good & if_qtl &
                         !duplicated(paste(chr, pos, sep = ":")) )
gwas$SNP_ID = paste(gwas$chr, gwas$pos, sep = ":")

if(nrow(gwas) == 0) stop("No individuals in the specified population are phenotyped for this GWAS trait!")


########## qtl regions with overlapped SNPs ##########
qtlColocReg_gwas = qtlColocReg[qtlColocReg$SNP_ID %in% gwas$SNP_ID, ]


########## coloc regions based on qtl ##########
gwas$chr = as.integer(gwas$chr)

gwasColocReg = qtlColocReg_gwas %>% select(c(Signal, Module, SNP_ID, Chr, Pos, Region)) %>% left_join(y = gwas, by = c("SNP_ID", "Chr"="chr", "Pos"="pos"))


# range of pvalues in a region for GWAS and qtl & remove the regions whose min GWAS p < 1e-5
gwasRegP = gwasColocReg %>% group_by(Region) %>% summarise(s = min(pval_EUR), l = max(pval_EUR))
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
