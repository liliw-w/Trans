########## Prepare coloc regions for gwas traits ##########
########## These traits are part of the 72 traits ##########
########## Focus on the 14 autoimmune diseases ##########
rm(list = ls()[!ls() %in% "qtlColocReg"])
library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

## qtl regions
file_qtlColocReg = "~/xuanyao_llw/coloc/qtlColocReg.txt.gz"

qtlColocReg = fread(file_qtlColocReg)


########## files and parameters ##########
gwas_pmid = 31604244
gwas_label = "MS"
p_included_thre = 1e-5

## create folder for each gwas trait
dir_gwas = file.path("/scratch/midway2/liliw1/coloc_MSigDB",
                     paste0("pmid", gwas_pmid, "_", gwas_label))
dir_gwas_data = file.path(dir_gwas, "data")
dir.create(dir_gwas, showWarnings = FALSE)
dir.create(dir_gwas_data, showWarnings = FALSE)

## gwas files
file_gwas = list.files(path = "/project2/xuanyao/llw/GWAS",
                       pattern = paste0("^pmid", gwas_pmid, ".*", gwas_label, ".*.gz$"),
                       full.names = TRUE)
file_gwas_trait_info = "/project2/xuanyao/llw/GWAS/72_traits_GWAS.csv"

## output files
file_qtlColocReg_gwas = file.path(dir_gwas_data, "qtlColocReg_gwas.txt.gz")
file_gwasColocReg = file.path(dir_gwas_data, "gwasColocReg.txt.gz")
file_gwasRegTruncPthre = file.path(dir_gwas_data, "gwasRegTruncPthre.txt")


########## read files ##########
## gwas meta files, info of 72 traits
gwas_trait_info = fread(file_gwas_trait_info)
gwas_trait_info[gwas_trait_info == ""] = NA

## gwas file
gwasCol = with(gwas_trait_info %>% filter(Label == gwas_label & PMID == gwas_pmid),
               c(snp_col_name, chr_col_name, pos_col_name,
                 af_col_name,
                 beta_col_name, se_col_name, p_col_name,
                 N_col_name, s_col_name))
gwasCol_newName = c("snp_rsID", "chr", "pos", "af", "beta", "se", "pval", "N", "s")
gwas = fread(file = file_gwas,
             select = gwasCol[!is.na(gwasCol)],
             col.names = gwasCol_newName[!is.na(gwasCol)])



########## gwas good SNPs & if overlapped with qtl SNPs & remove duplicated SNPs ##########
gwas$SNP_ID = with(gwas, paste(chr, pos, sep = ":"))

gwas$if_good = !is.na(gwas[["pval"]]) # & gwas$high_quality & !gwas[[paste0("low_confidence_", pop)]]
gwas$if_qtl = gwas$SNP_ID %in% qtlColocReg$SNP_ID
gwas$if_dup = duplicated(gwas$SNP_ID)
gwas = gwas %>% filter(if_good & if_qtl & !if_dup)

if(nrow(gwas) == 0) stop("No SNPs left for this GWAS trait after QC!")



########## Update qtl regions with overlapped SNPs ##########
ind_overlap = qtlColocReg$SNP_ID %in% gwas$SNP_ID

## Check if too few SNPs remained for coloc
cat("There were", nrow(qtlColocReg),
    "SNPs in the QTL coloc regions. \n",
    "After QC on GWAS SNPs and select overlapped SNPs, \n",
    sum(ind_overlap), "(", sum(ind_overlap)/nrow(qtlColocReg), ")",
    "SNPs are remained for coloc regions. \n")

qtlColocReg_gwas = qtlColocReg[qtlColocReg$SNP_ID %in% gwas$SNP_ID, ]


########## Construct gwas coloc regions based on qtl regions ##########
gwas$chr = as.numeric(gwas$chr)
gwasColocReg = qtlColocReg_gwas %>% select(c(Signal, Module, SNP_ID, Chr, Pos, Region, rsid)) %>% left_join(y = gwas, by = c("SNP_ID", "Chr"="chr", "Pos"="pos"))

## range of pvalues in a region for GWAS and qtl & remove the regions whose min GWAS p < 1e-5
gwasRegP = gwasColocReg %>% group_by(Region) %>% summarise(s = min(pval), l = max(pval))
qtlRegP = qtlColocReg_gwas %>% group_by(Region) %>% summarise(s = min(Pval), l = max(Pval))

## Regions whose min p < p_included_thre
gwasRegTruncPthre = gwasRegP$Region[gwasRegP$s < p_included_thre]


########## save results ##########
fwrite(qtlColocReg_gwas, file_qtlColocReg_gwas, quote = FALSE, sep = "\t")
fwrite(gwasColocReg, file_gwasColocReg, quote = FALSE, sep = "\t")
fwrite(data.table(gwasRegTruncPthre), file_gwasRegTruncPthre, quote = FALSE, sep = "\t", col.names = FALSE)

str(qtlColocReg_gwas)
str(gwasRegTruncPthre)
