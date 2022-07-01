############################
############## #SNP signals for each (module, chr) pair ##############
##########################################
rm(list = ls())
library(data.table)
library(tidyverse)


### read the indep_SNP file and construct df for plotting
df_indep_SNP <- fread("/scratch/midway2/liliw1/figures/fig3/fig3_data/num_indep_SNP.txt",
                      col.names = c('pair', 'num_indep_SNP'),
                      header = FALSE)
df_indep_SNP <- df_indep_SNP %>%
  separate(pair, c("module", "chr", NA), sep = '[.]') %>%
  separate(chr, c(NA, "chr"), sep = 'chr')


##########################################
############## plot Figure3B ##############
##########################################
### df for plotting
plot_df <- df_indep_SNP %>%
  separate(module, c(NA, "module_num"), "module", remove = FALSE)
# Chr order
plot_df$chr = factor(plot_df$chr, levels = 22:1)

# module order by #chr's on which there are signals & #indep signals
module_order <- plot_df %>%
  group_by(module, module_num) %>%
  summarise(n_M_chr = n(),
            n_M_indep_SNP = sum(num_indep_SNP) ) %>%
  ungroup() %>%
  arrange(desc(n_M_chr), desc(n_M_indep_SNP))
chr_order <- plot_df %>%
  group_by(chr) %>%
  summarise(n_M = n(),
            n_indep_SNP = sum(num_indep_SNP) ) %>%
  ungroup() %>%
  arrange(n_M, n_indep_SNP)
plot_df$module = factor(plot_df$module,
                        levels = module_order$module,
                        labels = module_order$module_num)
plot_df$chr = factor(plot_df$chr,
                     levels = chr_order$chr,
                     labels = chr_order$chr)


blue_seq <- c("#667ec6", "#4c68bd", "#3252b3", "#0028a1", "#002080", "#001860")


### plot
fig <- ggplot(plot_df, aes(x = chr, y = module)) +
  geom_tile(aes(fill = num_indep_SNP)) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(8, "Blues")[3:8]) +
  labs(x = "Chromosome", y = "Module", fill = "Independent Loci") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5, angle = -90, face = "bold"),
        legend.background = element_rect(color = "black",
                                         linetype = "dashed",
                                         size = 0.4),
        legend.key.size= unit(0.2, "cm"),
        axis.line = element_line(colour="black"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_line(size = 0.5),
        axis.text.x = element_text(angle = 90, colour = "black",
                                   vjust = 0, hjust = 1,
                                   size = 5),
        axis.text.y = element_text(colour = "black", size = 6.5),
        axis.title.x = element_text(vjust = -0.2, size=14),
        axis.title.y = element_text(angle=90,vjust =2, size=14) ) +
  guides(fill = guide_colourbar(title.position = "right")) +
  coord_flip()

fig

### save plot for further editing
saveRDS(fig, 'num_indep_SNP.rds')

### save plot
ggsave('num_indep_SNP.pdf',
       fig,
       width = 7, height = 3)
