rm(list = ls())
library(data.table)
library(tidyverse)
library(cowplot)


########## files and parameters, read files ##########
pvalThre = 'module_QTL_sig'
regionDis = 1e5
cisDis = 5e5
p_included_thre = 1e-5

file_module_QTL_signals = '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'
module_QTL_signals = fread(file_module_QTL_signals, col.names = c("signal", "p", "q"))
range(module_QTL_signals$p)
range(module_QTL_signals$q)

file_gene_meta = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
file_qtlColocReg = "~/xuanyao_llw/coloc/cis/data/qtlColocReg.txt.gz"

gene_meta = fread(file_gene_meta, header = TRUE)
qtlColocReg = fread(file_qtlColocReg)
qtlColocReg = qtlColocReg %>% filter(rsid != "")


### only remain the light green regions, i.e. regions whose lead-SNPs are trans-eQTL signals
res_coloc_leadSNP = qtlColocReg %>%
  group_by(Region) %>%
  summarise() %>%
  mutate(if_module_QTL = Region %in% module_QTL_signals$signal) %>%
  filter(if_module_QTL) %>%
  ungroup() %>%
  select(Region) %>%
  separate(col = Region, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)


########## loop over all traits ##########
#file_resColoc = paste0("/scratch/midway2/liliw1/coloc/data/pheno", gwasPhenocode_seq, ".resColoc.txt.gz")
#resColoc_all = rbindlist(lapply(file_resColoc, fread, header = TRUE), use.names = TRUE)

# only retain regions with module-QTL signals
#res_coloc_leadSNP = resColoc_all %>% filter(Pval <= pvalThre) %>%
#  select(c(Phenocode, trait, Region, Pval, nsnps, PP.H4.abf)) %>%
#  separate(col = Region, into = c("Module", "Chr", "Pos"), sep = ":", remove = FALSE, convert = TRUE)


# only consider "protein_coding" or "lincRNA" genes, & convert TSS
gene_meta = gene_meta %>% filter(Class %in% c("protein_coding", "lincRNA"))
strand_ind = gene_meta$Strand == "-"
gene_meta[strand_ind, "Start"] = gene_meta[strand_ind, "End"]


# TSS and distance of genes that are overlapped with regions and lead SNPs
coloc_cis_dis = max(regionDis/2, cisDis/2)

cis_overlap = rbindlist(apply(res_coloc_leadSNP, 1, function(x) {
  tmp_cis_overlap = gene_meta %>%
    filter(Chromosome == paste0("chr", as.numeric(x["Chr"]))) %>%
    mutate(dis = abs(Start-as.numeric(x["Pos"]))) %>%
    filter(dis < coloc_cis_dis)
  data.table("N_gene" = length(tmp_cis_overlap$GeneSymbol),
             "gene" = paste(tmp_cis_overlap$GeneSymbol, collapse = ";"),
             "dis" = paste(tmp_cis_overlap$dis, collapse = ";") )
}))


res_coloc_leadSNP = cbind(res_coloc_leadSNP, cis_overlap)


reg_gene = rbindlist(apply(res_coloc_leadSNP %>% filter(N_gene > 0), 1, function(x) {
  data.table("Region" = x["Region"],
             "gene" = unlist(strsplit(x["gene"], ";")) )
}))


#reg_gene[which(reg_gene$gene == 'ARHGEF3'), ]


file_cis_eqtl = "/project2/xuanyao/llw/DGN_PCO.lambda.01/nom_20phenoPC_allChr_eQTL.txt.gz"
cis_eqtl = fread(file_cis_eqtl)
colnames(cis_eqtl) = c("pid", "sid", "dist", "npval", "slope")


ind = cis_eqtl$pid %in% unique(reg_gene$gene)
cis_eqtl = cis_eqtl[ind, ]

#cis_eqtl[which(cis_eqtl$pid == 'ARHGEF3'), ]


gwasColocReg = reg_gene %>% left_join(cis_eqtl, by = c("gene" = "pid") ) %>% filter(!is.na(npval))
gwasColocReg = gwasColocReg %>%
  group_by(Region, gene) %>%
  summarise("minp" = min(npval)) %>%
  ungroup() %>%
  filter(minp < p_included_thre) %>%
  left_join(gwasColocReg, by = c("Region", "gene") )


gwasColocReg = gwasColocReg[gwasColocReg$sid %in% qtlColocReg$SNP_ID, ]
qtlColocReg_gwas = qtlColocReg[qtlColocReg$Region %in% gwasColocReg$Region & qtlColocReg$SNP_ID %in% gwasColocReg$sid, ]
gwasColocReg = gwasColocReg[gwasColocReg$sid %in% qtlColocReg_gwas$SNP_ID, ]

gwasColocReg = qtlColocReg_gwas %>%
  filter(!duplicated(SNP_ID)) %>%
  select(SNP_ID, MAF) %>%
  right_join(gwasColocReg, by = c("SNP_ID" = "sid") )


#gwasColocReg[which(gwasColocReg$gene == 'ARHGEF3'), ]

#qtlColocReg_gwas[qtlColocReg_gwas$Region =="module4:3:56901292", ]

fwrite(gwasColocReg, "gwasColocReg.txt.gz", quote = FALSE, sep = "\t")
fwrite(qtlColocReg_gwas, "qtlColocReg_gwas.txt.gz", quote = FALSE, sep = "\t")



########## Plot ##########
fig1 = ggplot(res_coloc_leadSNP, aes(x = N_gene)) +
  #facet_wrap(~trait) +
  geom_line(stat = "bin", binwidth = 1) +
  geom_point(stat = "bin", binwidth = 1, size = 0.7) +
  labs(x = "#cis genes", y = "#regions") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

res_coloc_dis_trait = res_coloc_leadSNP %>% group_by() %>% summarise(dis_trait = paste0(dis, collapse = ";"))
res_coloc_dis_trait = data.table("dis" = as.numeric(unlist(strsplit(res_coloc_dis_trait$dis_trait, ";"))) )

fig2 = ggplot(res_coloc_dis_trait, aes(x = dis)) +
  #facet_wrap(~trait) +
  geom_line(stat = "bin", binwidth = 5000) +
  geom_point(stat = "bin", binwidth = 5000, size = 0.7) +
  labs(x = "dis", y = "#regions") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

fig <- plot_grid(fig1, fig2, labels = c('A', "B"), nrow = 1)
ggsave("region.cis.gene.png", fig, path = "./plot/", width = 10, height = 5)


