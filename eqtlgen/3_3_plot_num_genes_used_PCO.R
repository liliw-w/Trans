############################################################
########### Plot the number of genes used for SNPs ###########
########### for each module ###########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
source("~/Trans/plot/theme_my_pub.R")

### I/O
file_num_gene_snp_used_all<- 'p/num_gene_snp_used_module_all.Sigma_nullz.txt'
num_gene_snp_used_all<- fread(file_num_gene_snp_used_all)


### convert column types
num_gene_snp_used_all <- num_gene_snp_used_all %>%
  type_convert(col_types = cols(module = col_double(),
                                module_size = col_double(),
                                SNP = col_character(),
                                module_size_eqtlgen = col_double(),
                                n_gene_trans = col_double(),
                                n_gene_trans_cross = col_double(),
                                meta = col_character()) )
module_seq <- num_gene_snp_used_all %>% distinct(module) %>% pull(module)

### plot for each module
for(module in module_seq){
  # extract for each module
  fig_dat <- num_gene_snp_used_all %>% filter(module == {{module}})
  
  # plot
  fig <- ggplot(fig_dat, aes(x = factor(SNPChr))) +
    geom_boxplot(aes(y = n_gene_trans_cross, color = "After cis+cross"), outlier.size = 0.3, show.legend = FALSE) +
    geom_jitter(aes(y = n_gene_trans_cross, color = "After cis+cross"), shape = 16, size = 0.3, alpha = 0.5) +
    #geom_quasirandom(aes(y = n_gene_trans_cross, color = "After cis+cross"), method = "tukey", shape = 16, size = 0.3, alpha = 0.5) +
    geom_line(aes(y = module_size, color = "Module size"), group = NA, linetype = "dashed") +
    #geom_line(aes(y = module_size_eqtlgen+1), group = NA, linetype = "dashed", color = "red") +
    geom_line(aes(y = n_gene_trans, color = "After cis"), group = NA, linetype = "dashed") +
    geom_point(aes(y = n_gene_trans, color = "After cis"), group = NA) +
    labs(x = "Chromosome", y = "Number of Genes", title = paste("Module", module) )
  
  # add themes
  fig +
    theme_my_pub(legend.position = "bottom") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed")) +
    scale_color_manual(name = "Gene type",
                       breaks = c("Module size", "After cis", "After cis+cross"),
                       values = c("Module size" = "black", "After cis" = "blue", "After cis+cross" = "#990000"),
                       guide = guide_legend() )
  
  # save
  ggsave(paste0("plot/num_gene_snp_used_M", module, ".pdf"), width = 6, height = 4)
  
  # verbose
  cat("Module", module, "is done.\n")
}

