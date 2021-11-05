rm(list = ls())
library(data.table)
library(tidyverse)


########## files and parameters ##########
file_gwas_bgz = "/scratch/midway2/liliw1/UKB_nealelab/continuous-30080-both_sexes-irnt.tsv.bgz"
file_snp_meta = "/scratch/midway2/liliw1/UKB_nealelab/full_variant_qc_metrics.txt.bgz"
file_qtlColocReg = "/scratch/midway2/liliw1/coloc/qtlColocReg.txt.gz"
file_qtlColocReg_gwas = "/scratch/midway2/liliw1/coloc/qtlColocReg_gwas.txt.gz"
file_gwasColocReg = "/scratch/midway2/liliw1/coloc/gwasColocReg.txt.gz"
file_gwasRegTruncPthre= "/scratch/midway2/liliw1/coloc/gwasRegTruncPthre.txt"

pop = "EUR"
regionDis = 1e5


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


########## qtl regions with overlapped SNPs ##########
qtlColocReg_gwas = qtlColocReg[qtlColocReg$SNP_ID %in% gwas$SNP_ID, ]


########## coloc regions based on qtl ##########
gwas$chr = as.integer(gwas$chr)

gwasColocReg = qtlColocReg_gwas %>% select(c(Signal, Module, SNP_ID, Chr, Pos, Region)) %>% left_join(y = gwas, by = c("SNP_ID", "Chr"="chr", "Pos"="pos"))


# range of pvalues in a region for GWAS and qtl & remove the regions whose min GWAS p < 1e-5
gwasRegP = gwasColocReg %>% group_by(Region) %>% summarise(s = min(pval_EUR), l = max(pval_EUR))
qtlRegP = qtlColocReg_gwas %>% group_by(Region) %>% summarise(s = min(Pval), l = max(Pval))

gwasRegTruncPthre = gwasRegP$Region[gwasRegP$s < 1e-5]


########## save results ##########
fwrite(qtlColocReg_gwas, file_qtlColocReg_gwas, quote = FALSE, sep = "\t")
fwrite(gwasColocReg, file_gwasColocReg, quote = FALSE, sep = "\t")
fwrite(gwasColocReg, file_gwasRegTruncPthre, quote = FALSE, sep = "\t", col.names = FALSE)


### Notes

# for trait 30080 Platelet count, there are 28,987,535 rows.
# For snp meta, there are 28,987,535 rows.

### only considering high-quality variants, `high_quality`
# out of 28987534 snps, 26691267 (~92%) are high_quality.

### retained variants with INFO scores > 0.8 (29,865,259 variants on the autosomes and X chromosome)
# For all snps in snp_meta file, the minimum info score is 0.8. So these snps are already filtered.

### low-confidence statistics removed, `low_confidence_{pop}`
