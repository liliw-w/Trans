############################################################
########### Plot #Indep null SNPs v.s. module size ###########
########### under a pre-spefified z-score threshold for being null ###########
############################################################
rm(list = ls())
library(tidyverse)
library(data.table)


### paras and I/O
thre_p_z <- 1e-4
file_null_SNP <- 'null_SNP/num_nullSNP.rds'
file_num_gene_used <- 'p/num_gene_snp_used_module_all.Sigma_nullz.txt'


# output
file_plot_nullSNP_module_size_line <- paste0('plot/plot_nullSNP_module_size_pthre', thre_p_z, '_ratio_chr.pdf')
file_ratio_module_chr <- 'null_SNP/ratio_module_chr.txt'


### read data
res_nullSNP <- readRDS(file_null_SNP)
num_gene_used <- fread(file_num_gene_used, header = TRUE)


fig_dat <- distinct(
  num_gene_used,
  module, SNPChr, module_size, module_size_eqtlgen, n_gene_trans
) %>%
  left_join(
    res_nullSNP %>%
      filter(thre_z == thre_p_z) %>%
      select(module, num_nullSNP_indep),
    by = "module"
  ) %>%
  mutate(
    prop_dim_module = num_nullSNP_indep/module_size_eqtlgen,
    prop_dim_module_chr = num_nullSNP_indep/n_gene_trans,
  )


ratio_thre <- c(5, 10, 50, 100, 150, Inf)
n_panel <- length(ratio_thre)
panel_title <- c(paste0(c(0, ratio_thre[-c(n_panel, n_panel-1)]),
                        "<Ratio<",
                        ratio_thre[-n_panel]),
                 paste0("Ratio>", ratio_thre[n_panel-1]))


fig_dat <- cbind(fig_dat,
                 "group_ratio" = apply(fig_dat, 1, function(x) which(as.numeric(x["prop_dim_module"]) <= ratio_thre)[1])
                 )
fig_dat$module <- factor(fig_dat$module, levels = 1:1000)
fig_dat$group_ratio <- factor(fig_dat$group_ratio,
                              levels = 1:n_panel,
                              labels = panel_title)



ggplot(fig_dat, aes(x = module)) +
  geom_point(aes(y = prop_dim_module),
             color = "#990000", size = 1) +
  geom_segment(aes(y = prop_dim_module, xend = module, yend = 0),
               color = "#990000", size = 0.3) +
  geom_boxplot(aes(y = prop_dim_module_chr, group = module),
               fill = NA, outlier.size = 0, show.legend = FALSE) +
  geom_jitter(aes(y = prop_dim_module_chr, group = module),
              shape = 16, size = 0.5, alpha = 0.5) +
  facet_wrap(~group_ratio, scales = "free") +
  labs(x = "Module", y = "Ratio -\n Null SNPs/Module size") +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.6, linetype = "solid"),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed", color = "#e5e5e5"),
        
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 90, size = 6),
        
        axis.ticks = element_blank(),
        
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "#f5f5f5", colour = "black", size = 0.6, linetype = "solid")
  )

ggsave(file_plot_nullSNP_module_size_line, height = 5, width = 8)


# save ratios for each module
fwrite(
  fig_dat,
  file_ratio_module_chr,
  quote = FALSE, sep = "\t"
)
