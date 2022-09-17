############################################################
########### Plot how many cross map genes are there for cis genes of each SNP ###########
############################################################
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggforce)

### I/O
file_num_cross_map_genes <- "num_cross_map_genes.rds"
file_cross_map_genes <- "cross_map_genes.rds"
file_plot_num_cross_map_genes <- 'plot/plot_num_cross_map_genes.pdf'
file_plot_num_cross_map_genes_order_by_chr <- 'plot/plot_num_cross_map_genes_order_by_chr.pdf'

### read files
num_cross_map_genes <- readRDS(file_num_cross_map_genes)
cross_map_genes <- readRDS(file_cross_map_genes)

### add SNP info to their numerical results of num of cross & cis genes
num_cross_map_genes <- cbind(cross_map_genes %>% select(1:4), num_cross_map_genes)



### 1. plot #cis-genes and #cross map genes for each SNP, ordered by num of cis-genes
fig_dat <- num_cross_map_genes %>%
  arrange(num_cis_gene) %>%
  mutate("num_cis_gene" = -num_cis_gene*100, "SNP" = row_number()) %>%
  pivot_longer(c(num_cis_gene, num_cross), names_to = "type", values_to = "num_gene")

fig <- ggplot(fig_dat, aes(x = SNP, y = num_gene, fill = type)) +
  geom_bar(stat = "identity", position = position_identity()) +
  geom_text(data = fig_dat %>% filter(type == "num_cis_gene" & between(SNP, 2000, 2020)),
            aes(y = num_gene, label = -num_gene/100),
            nudge_y = -200, size = 3) +
  facet_zoom(xy = between(SNP, 2000, 2020), zoom.size = 0.4) +
  labs(x = 'Chr/SNP', y = "Number of cross map genes with cis- genes \n Number of cis- genes", fill = " ") +
  coord_cartesian(ylim =  c(-2000, 6000))

fig +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        legend.position = "top") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Reds")[8],
                               RColorBrewer::brewer.pal(8, "Blues")[8]) )



ggsave(file_plot_num_cross_map_genes, height = 5, width = 10)



### 2. plot #cis-genes and #cross map genes for each SNP, ordered by chr and num of cis-genes
fig_dat <- num_cross_map_genes %>%
  arrange(SNPChr, num_cis_gene) %>%
  mutate("num_cis_gene" = -num_cis_gene*100, "SNP" = row_number()) %>%
  pivot_longer(c(num_cis_gene, num_cross), names_to = "type", values_to = "num_gene")
x_text_dat <- fig_dat %>% distinct(SNPChr, .keep_all = T)


fig <- ggplot(fig_dat, aes(x = SNP, y = num_gene, fill = type)) +
  geom_bar(stat = "identity", position = position_identity()) +
  geom_text(data = fig_dat %>% filter(type == "num_cis_gene" & between(SNP, 2000, 2020)),
            aes(y = num_gene, label = -num_gene/100),
            nudge_y = -200, size = 3) +
  facet_zoom(xy = between(SNP, 2000, 2020), zoom.size = 0.4) +
  labs(x = 'Chr/SNP', y = "Number of cross map genes with cis- genes \n Number of cis- genes", fill = " ") +
  coord_cartesian(ylim =  c(-2000, 6000))

fig +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        legend.position = "top") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Reds")[8],
                               RColorBrewer::brewer.pal(8, "Blues")[8]) ) +
  scale_x_continuous(breaks = pull(x_text_dat, SNP),
                     labels = pull(x_text_dat, SNPChr))


ggsave(file_plot_num_cross_map_genes_order_by_chr, height = 5, width = 10)
