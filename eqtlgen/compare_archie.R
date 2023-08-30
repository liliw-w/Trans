##############################################
########### Compare to archie results ###########
##############################################
# load packages -----
rm(list = ls())
library(tidyverse)
source('~/Trans/plot/theme_my_pub.R')

theme_set(
  theme_my_pub(legend.position = "right")
)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

# file_snp_all <- '/project2/xuanyao/llw/eQTLGen/eqtlgen_snps.txt'
# file_gene_all <- '/project2/xuanyao/llw/eQTLGen/eqtlgen_genes.txt'

file_signal_archie <- 'archie_results.xlsx'
file_signal_us <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/qtl_table.txt'


file_snp_used_us <- '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'
file_module_us <- '/project2/xuanyao/llw/DGN_no_filter_on_mappability/result/coexp.module.rds'
file_module_used_us <- '/project2/xuanyao/llw/eQTLGen/DGN_est_Sigma/postanalysis/module_use_ratio_50.txt'


## output -----
file_res <- 'compare_archie.txt'
file_fig <- 'fig_compare_archie.pdf'



# read files -----
# snp_all <- data.table::fread(file_snp_all)
# gene_all <- data.table::fread(file_gene_all)

signal_archie <- readxl::read_excel(file_signal_archie)
signal_us <- data.table::fread(file_signal_us)

snp_used_us <- data.table::fread(file_snp_used_us)
module_us <- readRDS(file_module_us)
module_used_us <- data.table::fread(file_module_used_us)



# organize data -----
# Of the 166 co-expression gene modules identified in DGN, 
# we used 129 modules with reliable correlation matrix approximations to ensure the trans-eQTL signals are well-controlled for inflation
# 
# We identified 8116 trans-eQTL SNP-gene coexpression module pairs, 
# corresponding to 2161 eQTLGen test SNPs and 122 gene modules (Figure 5B, Table S11); 
module_us <- enframe(module_us$moduleLabels, "gene", "m") %>%
  filter(m %in% module_used_us$module)

signal_archie <- group_by(signal_archie, Trait, ARCHIE_Component) %>%
  summarise(
    "Genes" = paste(Genes, collapse = ";"),
    "Variants" = paste(Variants[!is.na(Variants)], collapse = ";"),
    "Chrom_Position" = paste(Chrom_Position[!is.na(Chrom_Position)], collapse = ";")
  ) %>%
  ungroup()




# calculate number to compare signals -----
res <- cbind(
  select(signal_archie, Trait, ARCHIE_Component),
  
  apply(signal_archie, 1, function(x) {
    tmp_genes = strsplit(x['Genes'], ';') %>% unlist()
    tmp_var = strsplit(x['Variants'], ';') %>% unlist()
    tmp_chrpos = strsplit(x['Chrom_Position'], ';') %>% unlist()
    
    
    ## harmonize genes & variants -----
    tmp_m = module_us[match(tmp_genes, module_us$gene), ] %>% pull(m)
    if_m_used = !is.na(tmp_m)
    if_m_sig = tmp_m %in% signal_us$gene_module
    
    tmp_var = snp_used_us[match(tmp_var, snp_used_us$SNP), ] %>% pull(SNP)
    if_var_used = !is.na(tmp_var)
    if_var_sig = tmp_var %in% signal_us$rsid
    if_var_sig_m = tmp_var %in% 
      (filter(signal_us, gene_module %in% tmp_m) %>% pull(rsid))
    
    tibble(
      ## gene-wise -----
      'num_gene' = length(tmp_genes),
      
      # 1. Where are the selected genes distributed in gene modules?
      'sum_m_used' = sum(if_m_used),
      
      
      # 2. If they are in trans target modules?
      ## if selected genes are replicated in trans target gene modules
      'sum_m_sig' = sum(if_m_sig),
      
      
      
      
      ## var-wise -----
      'num_var' = length(tmp_var),
      
      # 3. If selected variants are included in analysis?
      'sum_var_used' = sum(if_var_used),
      
      # 4. If selected variants are replicated in trans signals?
      'sum_var_sig' = sum(if_var_sig),
      
      
      # 5. If selected variants are signals for target module selected genes in
      'sum_var_sig_m' = sum(if_var_sig_m),
      
      ## misc -----
      'if_m_used' = paste(if_m_used, collapse = ";"),
      'if_m_sig' = paste(if_m_sig, collapse = ";"),
      'if_var_sig' = paste(if_var_sig, collapse = ";"),
      'if_var_sig_m' = paste(if_var_sig_m, collapse = ";"),
      'tmp_m' = paste(tmp_m, collapse = ";")
    )
  }) %>%
    bind_rows()
)



# visualize -----
ggpubr::ggarrange(
  ggplot(res, aes(x = str_glue('{Trait}_{ARCHIE_Component}'))) +
    geom_col(aes(y = num_gene, fill = "num_all")) +
    geom_col(aes(y = sum_m_used, fill = "sum_used")) +
    geom_col(aes(y = sum_m_sig, fill = "sum_sig")) +
    labs(x = NULL, y = "Number of selected genes", fill = 'Label') +
    scale_fill_manual(
      values = c("num_all" = "#C6DBEF", "sum_used" = "#6BAED6", "sum_sig" = "#2171B5"),
      breaks = c("num_all", "sum_used", "sum_sig"),
      labels = c("num_all" = "Selected", "sum_used" = "Included", "sum_sig" = "Signal")
    ),
  
  ggplot(res, aes(x = str_glue('{Trait}_{ARCHIE_Component}'))) +
    geom_col(aes(y = num_var, fill = "num_all")) +
    geom_col(aes(y = sum_var_used, fill = "sum_used")) +
    geom_col(aes(y = sum_var_sig, fill = "sum_sig")) +
    geom_errorbar(aes(y = sum_var_sig_m, ymin = sum_var_sig_m, ymax = sum_var_sig_m), color = "maroon") +
    labs(x = NULL, y = "Number of selected variants", fill = 'Label') +
    scale_fill_manual(
      values = c("num_all" = "#C6DBEF", "sum_used" = "#6BAED6", "sum_sig" = "#2171B5"),
      breaks = c("num_all", "sum_used", "sum_sig"),
      labels = c("num_all" = "Selected", "sum_used" = "Included", "sum_sig" = "Signal")
    ),
  
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  common.legend = TRUE, legend = "right"
)



# print out key message or write out -----
data.table::fwrite(res, file_res, sep = "\t")

ggsave(
  filename = file_fig,
  width = 6, height = 4
)

