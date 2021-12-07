###################################################
#######Check if lead SNPs of coloc regions are eQTLGen signals###########
###################################################
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)


########## files and parameters, read files ##########
pp4Thre = 0.75
pvalThre = 1e-8
nsnpsThre = 5

gwasPhenocode_seq = c(30080, 30090, 30100, 30110, 30010, 30020, 30030, 30040, 30050, 30060, 30070, 30270, 30240, 30250, 30260, 30280, 30290, 30300, 30000, 30120, 30130, 30140, 30150, 30160, 30180, 30190, 30200, 30210, 30220)


########## loop over all traits, for coloc regions ##########
file_resColoc = paste0("/scratch/midway2/liliw1/coloc/data/pheno", gwasPhenocode_seq, ".resColoc.txt.gz")
resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

## coloc regions
res_coloc_leadSNP = resColoc_all %>% filter(Pval <= pvalThre & PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) %>%
  select(c(Phenocode, trait, Region, Pval, nsnps, PP.H4.abf)) %>%
  separate(col = Region, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE) %>%
  mutate("SNP_id" = paste(Chr, Pos, sep = ":"))


file_trait_info = "/scratch/midway2/liliw1/coloc/ukbb_blood_traits.csv"
trait_info = fread(file_trait_info, sep = ",", header = TRUE)

module = 4

res_coloc_leadSNP = res_coloc_leadSNP %>%
  filter(Module == paste0("module",module)) %>%
  left_join(trait_info, by = c(Phenocode = 'GWAS ID', trait = 'Trait Abbreviation') ) %>%
  arrange(SNP_id) %>%
  select(c(Phenocode, trait, 'GWAS Group', 'GWAS Trait', Region))

fwrite(res_coloc_leadSNP, paste0("module",module,".coloc.traits.txt"), quote = FALSE, sep = "\t")


df1 = as.data.frame(table(res_coloc_leadSNP$Region, dnn = list("Region")), responseName = "traits")
df1 = df1 %>% separate(col = Region, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)

fwrite(df1, paste0("module",module,".coloc.traits.snp.summary.txt"), quote = FALSE, sep = "\t")


file_gene_meta = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
gene_meta = fread(file_gene_meta, header = TRUE)
#gene_meta = gene_meta %>% filter(Class %in% c("protein_coding", "lincRNA"))
strand_ind = gene_meta$Strand == "-"
gene_meta[strand_ind, "Start"] = gene_meta[strand_ind, "End"]


for(i in 1:nrow(df1)){
  x = df1[i, ]
  tmp = gene_meta[gene_meta$Chromosome == paste0("chr", x['Chr']), ] %>% mutate("dis" = as.numeric(x["Pos"])-Start) %>% filter(abs(dis) < 5e+5/2) %>% arrange(abs(dis))
  print(x)
  print(tmp)
}



########## eQTLGen signals ##########
## eQTLGen files
file_eqtlgen_trans = "/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz"
file_eqtlgen_cis = "/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
file_snp_trans = "/project2/xuanyao/data/eQTLGen/meta.snp.txt.gz"
file_snp_cis = "/project2/xuanyao/data/eQTLGen/meta.snp.cis.txt.gz"

## read files
eqtlgen_trans = fread(file_eqtlgen_trans)
eqtlgen_cis = fread(file_eqtlgen_cis)
snp_trans = fread(file_snp_trans)
snp_cis = fread(file_snp_cis)

## mutate SNP_id
eqtlgen_trans$SNP_id = paste(eqtlgen_trans$SNPChr, eqtlgen_trans$SNPPos, sep = ":")
snp_trans$SNP_id = paste(snp_trans$SNPChr, snp_trans$SNPPos, sep = ":")
eqtlgen_cis$SNP_id = paste(eqtlgen_cis$SNPChr, eqtlgen_cis$SNPPos, sep = ":")
snp_cis$SNP_id = paste(snp_cis$SNPChr, snp_cis$SNPPos, sep = ":")


########## check coloc regions that are overlapped with eQTLGen signals ##########
N_region = nrow(res_coloc_leadSNP)
N_uniq_region = length(unique(res_coloc_leadSNP$Region))
N_uniq_leadSNP = length(unique(res_coloc_leadSNP$SNP_id))
N_uniq_leadSNP_trans_snp = sum(unique(res_coloc_leadSNP$SNP_id) %in% snp_trans$SNP_id)
N_uniq_leadSNP_trans_sig = sum(unique(res_coloc_leadSNP$SNP_id) %in% eqtlgen_trans$SNP_id)
N_uniq_leadSNP_cis_snp = sum(unique(res_coloc_leadSNP$SNP_id) %in% snp_cis$SNP_id)
N_uniq_leadSNP_cis_sig = sum(unique(res_coloc_leadSNP$SNP_id) %in% eqtlgen_cis$SNP_id)
N_uniq_leadSNP_eqtlgen_snp = sum(unique(res_coloc_leadSNP$SNP_id) %in% c(snp_trans$SNP_id, snp_cis$SNP_id) )
N_uniq_leadSNP_eqtlgen_sig = sum(unique(res_coloc_leadSNP$SNP_id) %in% c(eqtlgen_trans$SNP_id, eqtlgen_cis$SNP_id) )

## trans-
cat("#all regions for all traits:", N_region, ". \n",
    "#unique regions:", N_uniq_region, ". \n",
    "#unique lead SNPs:", N_uniq_leadSNP, ". \n",
    "Out of", N_uniq_leadSNP, "unique lead SNPs,", N_uniq_leadSNP_trans_snp, "are included in eQTLGen trans- SNPs. \n",
    "Out of", N_uniq_leadSNP_trans_snp, "SNPs,", N_uniq_leadSNP_trans_sig, "are eQTLGen trans- signals. \n")

## cis-
cat("#all regions for all traits:", N_region, ". \n",
    "#unique regions:", N_uniq_region, ". \n",
    "#unique lead SNPs:", N_uniq_leadSNP, ". \n",
    "Out of", N_uniq_leadSNP, "unique lead SNPs,", N_uniq_leadSNP_cis_snp, "are included in eQTLGen cis- SNPs. \n",
    "Out of", N_uniq_leadSNP_cis_snp, "SNPs,", N_uniq_leadSNP_cis_sig, "are eQTLGen cis- signals. \n")

## trans- or cis-
cat("Out of", N_uniq_leadSNP, "unique lead SNPs,", N_uniq_leadSNP_eqtlgen_snp, "are included in eQTLGen SNPs. \n",
    "Out of", N_uniq_leadSNP_eqtlgen_snp, "SNPs,", N_uniq_leadSNP_eqtlgen_sig, "are eQTLGen cis- or trans- signals. \n")

