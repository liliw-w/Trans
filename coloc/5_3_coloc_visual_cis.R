########## visualization of colocalization events ##########
########## by R package locuscomparer ##########
rm(list = ls())
library(data.table)
library(tidyverse)
library(locuscomparer)
library(cowplot)


## output files
file_qtlColocReg_gwas = "/project2/xuanyao/llw/coloc/cis/cis_e/data/qtlColocReg_gwas.txt.gz"
file_gwasColocReg = "/project2/xuanyao/llw/coloc/cis/cis_e/data/gwasColocReg.txt.gz"
file_resColoc_all = "/project2/xuanyao/llw/coloc/cis/cis_e/data/resColoc.txt.gz"

pp4Thre = 0.75
nsnpsThre = 5
popColoc = "EUR"
refGenom = 'hg19'


########## files and parameter, read files ##########
qtlColocReg_gwas = fread(file_qtlColocReg_gwas)
gwasColocReg = fread(file_gwasColocReg)
resColoc_all = fread(file_resColoc_all)

dir.create(file.path("colocVisData"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("plot"), recursive = TRUE, showWarnings = FALSE)


########## only retain regions with high pp4 ##########
#removeReg = c("module4:12:8142264")
#resColoc = resColoc %>% filter(PP.H4.abf > pp4Thre)
resColoc_all = filter(resColoc_all, PP.H4.abf > pp4Thre & nsnps >= nsnpsThre)


## remove SNPs without rsid, as locuscomparer requires rsid
qtlColocReg_gwas = qtlColocReg_gwas[qtlColocReg_gwas$rsid != "", ]
gwasColocReg$rsid <- qtlColocReg_gwas[match(gwasColocReg$SNP_ID, qtlColocReg_gwas$SNP_ID), ]$rsid
gwasColocReg = gwasColocReg[gwasColocReg$rsid != "", ]


gwasPhenocode_seq = unique(resColoc_all$gene)
gwas_trait_seq = unique(resColoc_all$gene)
n_gwas = length(unique(resColoc_all$gene))


# i <- 57, 196, 264, 272, 396:398 ??? "There must be one and only one chromosome."
for(i in 1:n_gwas){
  gwasPhenocode = gwasPhenocode_seq[i]
  gwas_trait = gwas_trait_seq[i]
  


  resColoc <- resColoc_all %>% filter(gene == gwasPhenocode)
  
  ########## Visualize each region & save ##########
  nRegion = length(resColoc$Region)
  
  for(reg in resColoc$Region){
    M = strsplit(reg, ":") %>% unlist() %>% .[1] %>% strsplit(., "module") %>% unlist() %>% .[2] %>% paste("Module", .)
    
    tmpresColoc = resColoc %>% filter(Region == reg)
    tmpqtlColocReg_gwas = qtlColocReg_gwas %>% filter(Region == reg)
    tmpgwasColocReg = gwasColocReg %>%
      filter(gene == gwasPhenocode &
               Region == reg &
               SNP_ID %in% tmpqtlColocReg_gwas$SNP_ID)
    
    file_qtlVis = paste0("colocVisData/", reg, ".", gwasPhenocode, ".qtl.reg.coloc.txt")
    file_gwasVis = paste0("colocVisData/", reg, ".", gwasPhenocode, ".gwas.reg.coloc.txt")
    
    fwrite(tmpqtlColocReg_gwas %>% select(rsid, Pval) %>% rename(pval = Pval),
           file_qtlVis, quote = FALSE, sep = "\t")
    fwrite(tmpgwasColocReg %>% select(rsid, npval) %>% rename(pval = npval),
           file_gwasVis, quote = FALSE, sep = "\t")
    
    figTitle = paste0("Region:", tmpresColoc$Region,
                      "; Gene:", gwasPhenocode,
                      "\n#SNPs:", tmpresColoc$nsnps,
                      "; PP.H4.abf:", format(tmpresColoc$PP.H4.abf, digits = 4))
    basic_plt <- locuscompare(in_fn1 = file_gwasVis, in_fn2 = file_qtlVis,
                        title1 = gwasPhenocode, title2 = M,
                        population = popColoc, genome = refGenom)
    
    fig <- basic_plt +
      labs(title = figTitle) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10)
      )
    
    saveRDS(fig, paste0("plot/", reg, ".", gwasPhenocode, ".reg.coloc.rds"))
    ggsave(paste0("plot/", reg, ".", gwasPhenocode, ".reg.coloc.pdf"), fig, width = 4, height = 2.5)
    
    ### in progress
    cat(grep(reg, resColoc$Region, fixed = TRUE),
        "-th region:", reg, "(out of", nRegion, "regions)", "is done! \n",
        "Gene: ", gwasPhenocode, '\n')
  }
}


### Look at how many regions have coloc cis- genes; coloc to how many cis- genes; what are these genes
nsnpsThre <- 5
pvalThre <- "module_QTL"

res_coloc_reg_prop = resColoc_all %>%
  group_by(Region) %>%
  summarise(nGene = n(),
            nGeneColoc = sum(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) ) %>%
  ungroup()
res_coloc_reg_prop <- resColoc_all %>%
  filter(PP.H4.abf > pp4Thre & nsnps >= nsnpsThre) %>%
  group_by(Region) %>%
  summarise(gene = paste(gene, collapse = ";") ) %>%
  right_join(res_coloc_reg_prop, by = "Region") %>%
  arrange(desc(nGeneColoc), desc(nGene) )

res_coloc_reg_prop = res_coloc_reg_prop %>%
  mutate(propColoc = nGeneColoc/nGene )


fwrite(res_coloc_reg_prop,  paste0("data/coloc_region_prop_pvalThre-", pvalThre, ".txt"),
       quote = FALSE, sep = "\t")

