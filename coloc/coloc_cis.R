########## Run coloc on regions for qtl and gwas traits ##########
########## These traits are part of the 72 traits ##########
########## Focus on the 14 autoimmune diseases ##########
rm(list = ls())
library(data.table)
library(tidyverse)
library(coloc)

#setwd('/scratch/midway2/liliw1/coloc_cis_ARHGEF3')


## region files
file_qtlColocReg =  "~/xuanyao_llw/coloc/cis/qtlColocReg.txt.gz"
file_qtlColocReg_gwas = file.path("/scratch/midway2/liliw1/coloc_cis_eQTL/data/qtlColocReg_gwas.txt.gz")
file_gwasColocReg = file.path("/scratch/midway2/liliw1/coloc_cis_eQTL/data/gwasColocReg.txt.gz")


## output
file_resColoc = file.path("data/resColoc.txt.gz")
file_resColocSNP = file.path("data/resColocSNP.txt.gz")


qtlColocReg <- fread(file_qtlColocReg, header = TRUE)
qtlColocReg_gwas <- fread(file_qtlColocReg_gwas)
gwasColocReg <- fread(file_gwasColocReg)


gene_sep <- gwasColocReg %>% distinct(gene) %>% pull()
gwas_pmid_seq <- gene_sep
gwas_label_seq <- gene_sep


qtlN = 913
gwasN <- 913
qtlType = "quant"
gwasType = "quant"
pp4Thre = 0.75
csThre = 0.95

gwasColocReg_all_gene <- gwasColocReg
resColoc = NULL
resColocSNP = NULL
for(k in 1:length(gwas_pmid_seq)){
  gwas_pmid = gwas_pmid_seq[k]
  gwas_label = gwas_label_seq[k]
  
  gwasColocReg <- gwasColocReg_all_gene %>% filter(gene == gwas_pmid)
  
  gwasRegTruncPthre <- gwasColocReg %>% distinct(Region)
  
  
  ########## prepare coloc files and run coloc ##########
  nRegion = length(gwasRegTruncPthre$Region)
  for(reg in gwasRegTruncPthre$Region){
    ### extract the region for qtl and gwas
    tmpqtlColocReg_gwas = qtlColocReg_gwas %>% filter(Region == reg)
    tmpgwasColocReg = gwasColocReg %>% filter(Region == reg)
    
    tmpgwasColocReg <- tmpgwasColocReg %>% filter(!duplicated(SNP_ID)) %>% filter(SNP_ID %in% tmpqtlColocReg_gwas$SNP_ID)
    tmpqtlColocReg_gwas <- tmpqtlColocReg_gwas %>% filter(SNP_ID %in% tmpgwasColocReg$SNP_ID)
    
    ### construct coloc dataset D1 & D2
    D1 = list("pvalues" = tmpqtlColocReg_gwas$Pval,
              "N" = qtlN,
              "MAF" = tmpqtlColocReg_gwas$MAF,
              "type" = qtlType,
              "snp" = tmpqtlColocReg_gwas$SNP_ID)
    D2 = list("pvalues" = tmpgwasColocReg$npval,
              "N" = gwasN,
              "MAF" = tmpgwasColocReg$MAF,
              "beta" = tmpgwasColocReg$slope,
              #"varbeta" = (tmpgwasColocReg[["se"]])^2,
              "type" = gwasType,
              #"s" = if(gwasType=="cc" & !is.na(gwasN) ) n_cases/gwasN else if(gwasType=="cc" & is.na(gwasN) ) tmpgwasColocReg[["s"]] else NULL,
              "snp" = tmpgwasColocReg$SNP_ID)
    
    ### do coloc
    coloc_res = coloc.abf(D1, D2)
    
    ### follow-up: credible set
    if(coloc_res$summary["PP.H4.abf"] > pp4Thre){
      o <- order(coloc_res$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(coloc_res$results$SNP.PP.H4[o])
      w <- which(cs > csThre)[1]
      resColocSNP = rbind(resColocSNP, data.table("Region" = reg,
                                                  "SNP_ID" = as.character(coloc_res$results[o,][1:w,]$snp),
                                                  "SNP.PP.H4" = coloc_res$results[o,][1:w,]$SNP.PP.H4))
    }
    
    ### organize result
    resColoc = rbind(resColoc, data.table("Region" = reg, "gwas_label" = gwas_label, t(coloc_res$summary)) )
    
    ### in progress
    cat(grep(reg, gwasRegTruncPthre$Region, fixed = TRUE),
        "-th region:", reg, "(out of", nRegion, "regions)", "is done!", "\n")
  }
  
  cat("Trait", gwas_label, ". \n")
  
}


########## Add p-value and trait info to coloc regions ##########
resColoc = qtlColocReg %>%
  select(c(Signal, Pval)) %>%
  right_join(y = resColoc, by = c("Signal" = "Region")) %>%
  rename("Region" = "Signal")
resColoc <- resColoc %>%
  mutate(gene = gwas_label) %>%
  rename("Phenocode" = "gwas_label")


########## save results##########
resColoc = resColoc %>% arrange(desc(PP.H4.abf))
fwrite(resColoc, file_resColoc, quote = FALSE, sep = "\t")
#if(!is.null(resColocSNP)) fwrite(resColocSNP, file_resColocSNP, quote = FALSE, sep = "\t")


str(resColoc)
resColoc %>% distinct(Region)
resColoc %>% group_by(Region) %>%
  summarise(H0 = max(PP.H0.abf),
            H1 = max(PP.H1.abf),
            H2 = max(PP.H2.abf),
            H3 = max(PP.H3.abf),
            H4 = max(PP.H4.abf) ) %>%
  arrange(desc(H4)) %>%
  filter(
    H4  > pp4Thre
  )





df_fig <- resColoc %>% group_by(Region) %>% summarise(H0 = max(PP.H0.abf),
                                                      H1 = max(PP.H1.abf),
                                                      H2 = max(PP.H2.abf),
                                                      H3 = max(PP.H3.abf),
                                                      H4 = max(PP.H4.abf) ) %>%
  pivot_longer(c('H0', 'H1', 'H2', 'H3', 'H4'), names_to = "type", values_to = "H")

df_fig2 <- resColoc %>%
  pivot_longer(c('PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf'), names_to = "type", values_to = "H")


library(ggbeeswarm)
library(RColorBrewer)

custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(df_fig2, mapping = aes(type, H, color = type)) +
  geom_quasirandom(varwidth = TRUE, size = 1, alpha = 0.6) +
  scale_color_manual(values = cbp1) +
  theme_classic()


# 255
# 212


