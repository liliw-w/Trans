############################################################
########### Plot #Indep null SNPs (and #Uniq null SNPs) for each module ###########
########### under various z-score threshold for being null ###########
############################################################

rm(list = ls())
library(ggplot2)

## I/O
file_null_SNP <- 'null_SNP/num_nullSNP.rds'
file_plot_null_SNP <- 'null_SNP/num_nullSNP.pdf'
source('~/Trans/plot/theme_my_pub.R')

res_nullSNP <- readRDS(file_null_SNP)


## bar plot
base_fig <- ggplot(data = res_nullSNP, aes(x = factor(module)) ) +
  geom_bar(aes(y = -module_size), stat="identity") +
  geom_point(aes(y = num_nullSNP_uniq, color = factor(thre_z), shape = factor("Uniq null SNP")), alpha = 0.9) +
  geom_point(aes(y = num_nullSNP_indep, color = factor(thre_z), shape = factor("Indep null SNP")) ) +
  geom_text(data = res_nullSNP[!duplicated(res_nullSNP$module), ],
            aes(y = -module_size*10-200, label = module_size),
            size = 2) +
  scale_shape_manual(values = c(19, 3)) +
  labs(x = "Module", y = "Number of null SNP", color="Threshold of z", shape = "SNP Types")

## adjust themes
fig <- base_fig +
  theme_my_pub(axis.text.x.angle = 45) +
  theme(axis.text.x = element_text(size = 7),
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Dark2")


ggsave(file_plot_null_SNP, fig, height = 8, width = 12)

