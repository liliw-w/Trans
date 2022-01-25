rm(list = ls())
library(data.table)
library(tidyverse)
library(readxl)
library(ggplot2)

file_sup_tab8 = "/scratch/midway2/liliw1/Franke_2021.xlsx"
file_qtlColocReg = "/project2/xuanyao/llw/coloc/qtlColocReg.txt.gz"

trait_of_interest = "Systemic lupus erythematosus"
min_NrOfUnlinkedLoci = 1


sup_tab8 = read_xlsx(file_sup_tab8, skip = 3)
qtlColocReg = fread(file_qtlColocReg)

trait_all = unique(sup_tab8$Trait)

res_reg_Franke = NULL
for (trait_of_interest in trait_all) {
  df_trait_of_interest = sup_tab8 %>%
    filter(Trait == trait_of_interest & NrOfUnlinkedLoci >= min_NrOfUnlinkedLoci)
  
  tmp_snp_1 = unlist(strsplit(df_trait_of_interest$LociWithInformationOnIndividualTransEQTLs, ",", fixed = TRUE))
  tmp_snp_2 = unlist(strsplit(tmp_snp_1, "|", fixed = TRUE))
  tmp_snp_3 = unlist(strsplit(tmp_snp_2, "[", fixed = TRUE))
  signal_of_interest = unique(tmp_snp_3[grepl("^rs", tmp_snp_3)])
  
  qtlColocReg_signal = qtlColocReg %>% filter(rsid %in% signal_of_interest)
  tmp_res_reg_Franke = qtlColocReg_signal %>%
    select(c("rsid", "SNP_ID", "Region")) %>%
    mutate("trait_Franke" = trait_of_interest)
  res_reg_Franke = rbind(res_reg_Franke, tmp_res_reg_Franke)
  
}

file_resColoc_all = list.files("/project2/xuanyao/llw/coloc",
                               "resColoc.txt.gz",
                               full.names = TRUE, recursive = TRUE)
file_resColoc_all = file_resColoc_all[!grepl("cis", file_resColoc_all)]
res_Franke_coloc = res_reg_Franke
file_coloc_fig = NULL
for(file_resColoc in file_resColoc_all){
  resColoc = fread(file_resColoc)
  tmp_col = as.character(resColoc$Phenocode[1])
  res_Franke_coloc = resColoc %>%
    select(c("Region", "PP.H4.abf")) %>%
    right_join(y = res_Franke_coloc, by = c("Region")) %>%
    rename(!!paste0("PP4_", tmp_col) := PP.H4.abf)
}



saveRDS(res_Franke_coloc, "res_Franke_coloc.rds")

res_Franke_coloc = readRDS("res_Franke_coloc.rds")

for(trait_of_interest in trait_all){
  tmp = res_Franke_coloc %>%
    filter(trait_Franke == trait_of_interest) %>%
    pivot_longer(starts_with("PP4"), names_to = "trait_coloc", values_to = "PP4")
  tmp[is.na(tmp)] = 0
  tmp$Region = paste0(tmp$Region, "\n(", tmp$rsid, ",", tmp$SNP_ID, ")")
  
  fig_tile <- ggplot(tmp, aes(trait_coloc, Region)) +
    geom_tile(aes(fill = PP4)) +
    scale_fill_gradient(low="white", high="grey2") +
    labs(fill = "PP4") +
    geom_point(data = tmp[tmp$PP4 >= 0.75, ], shape = 8, color = "red", size = 0.5) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top") +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 5),
          axis.text.y = element_text(size = 5),
          axis.title=element_text(size=5))
  ggsave(paste0(trait_of_interest, "_tile.png"), fig_tile, path = "~/scratch/")
  
}

ggsave(paste0(trait_of_interest, "_tile.png"), fig_tile, path = "~/scratch/", height = 2.5, width = 7)

tmp = res_Franke_coloc %>%
  #filter(trait_Franke == trait_of_interest) %>%
  pivot_longer(starts_with("PP4"), names_to = "trait_coloc", values_to = "PP4")
tmp[is.na(tmp)] = 0
tmp$Region = paste0(tmp$Region, "\n(", tmp$SNP_ID, ")")

fig_tile <- ggplot(tmp, aes(trait_coloc, Region)) +
  geom_tile(aes(fill = PP4)) +
  scale_fill_gradient(low="white", high="grey2") +
  labs(fill = "PP4") +
  geom_point(data = tmp[tmp$PP4 >= 0.75, ], shape = 8, color = "red", size = 0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title=element_text(size=5))
ggsave(paste0("tmp_tile.png"), fig_tile, path = "~/scratch/", height = 30)



a = sup_tab8 %>% filter(Trait == "Systemic lupus erythematosus")
b = unique(a[grep("rs4917014", a$LociWithInformationOnIndividualTransEQTLs), ]$LociConvergeOn)

meta = fread("~/scratch/DGN_no_filter_on_mappability/result/gene.meta.txt")
meta$gene_ensg = sapply(meta$GeneNameConv, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
sum(b %in% meta$gene_ensg)
d = meta[meta$gene_ensg %in% b, ]$gene

coexp = readRDS("~/scratch/DGN_no_filter_on_mappability/result/coexp.module.rds")$moduleLabels
coexp[d]



