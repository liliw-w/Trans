############## plot Figure3B ##############
############## #SNP signals for each (module, chr) pair ##############
##########################################
rm(list = ls())
library(data.table)
library(tidyverse)

file_signal <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt'

signal <- fread(file_signal, col.names = c("module_SNP", "p", "q") )


### write unique SNP signals for each (module, chr) pair to file for follow up plink
signal %>%
  separate(module_SNP, c("module", "chr", "pos"), ":", remove = FALSE, convert = TRUE) %>%
  unite("SNP_ID", chr, pos, sep = ":", remove = FALSE) %>%
  select(module, chr, SNP_ID) %>%
  group_by(module, chr) %>%
  group_walk(~ fwrite(.x,
                      file.path("/scratch/midway2/liliw1/figures",
                                paste0(.y$module, ".chr", .y$chr, ".txt")),
                      quote = FALSE, col.names = FALSE) )

### Use plink to obtain independent SNP signals for each (module, chr) pair
df_indep_SNP = NULL
file_uniq_list <- list.files(path = "/scratch/midway2/liliw1/figures/fig3/fig3_data/",
                             pattern = "^module.*chr.*.txt$")
for(file_uniq in file_uniq_list){
  file_indep = paste0("indep.", file_uniq)
  command = paste0("bash /home/liliw1/xuanyao_llw/eQTLGen_est_Sigma/indep.SNP.sh ", file_uniq, " ", file_indep)
  system(command)
  df_indep_SNP = rbind(df_indep_SNP, c(
    file_uniq,
    as.numeric(system(paste0("cat ", file_indep, " | wc -l"), intern = TRUE))
  ))
}
fwrite(df_indep_SNP, "num_indep_SNP.txt", quote = FALSE, col.names = FALSE, sep = '\t')


### read the indep_SNP file and construct df for plotting
df_indep_SNP <- fread("./fig3/fig3_data/num_indep_SNP.txt",
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
plot_df$module = factor(plot_df$module,
                        levels = module_order$module,
                        labels = module_order$module_num)


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
        legend.title = element_text(size = 5, angle = 90, face = "bold"),
        legend.background = element_rect(color = "black",
                                         linetype = "dotted",
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
        axis.title.x = element_text(vjust = -0.2, size=9),
        axis.title.y = element_text(angle=90,vjust =2, size=9) ) +
  guides(fill = guide_colourbar(title.position = "left")) +
  coord_flip()


### save plot for further editing
saveRDS(fig, 'fig3b_num_indep_SNP.rds')

### save plot
ggsave('fig3b_num_indep_SNP.png',
       fig,
       width = 7, height = 3)
