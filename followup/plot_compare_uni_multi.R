### Plot the uni SNPs and those are also multi SNPs
rm(list = ls())
library(data.table)
library(tidyverse)

file_battle_QTL <- '/project2/xuanyao/llw/DGN_PCO.lambda.01/DGN_transeQTL_battle2014.txt'
file_module_QTL_signals <- '~/xuanyao_llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr5.txt'
file_coexp_module <- '~/xuanyao_llw/DGN_no_filter_on_mappability/result/coexp.module.rds'


battle_QTL <- fread(file_battle_QTL, header = TRUE)
module_QTL_signals <- fread(file_module_QTL_signals, header = FALSE, col.names = c("signal", "p", "q"))
coexp_module <- readRDS(file_coexp_module)$moduleLabels

coexp_module <- tibble(module = coexp_module, gene = names(coexp_module))
dgn_gene <- coexp_module$gene

module_QTL_signals <- module_QTL_signals %>%
  separate(signal, c('module', 'chr', 'pos'), ":", remove = FALSE) %>%
  unite(col = 'SNP_ID', c(chr, pos), sep = ":", remove = FALSE)

# how many unique genes
battle_QTL %>% distinct(GENE_NAME)

# how many unique SNPs
battle_QTL %>% distinct(SNP_ID)

# battle genes also considered in multi test pipeline
battle_QTL <- filter(battle_QTL, GENE_NAME %in% dgn_gene)

# how many unique genes & SNPs
battle_QTL %>% distinct(GENE_NAME)
battle_QTL %>% distinct(snp)


# how many unique SNPs, and their target genes
battle_QTL_snp <- battle_QTL %>% group_by(snp) %>% summarise(n_gene = n() )
battle_QTL_snp <- battle_QTL_snp %>% mutate(if_rep = snp %in% unique(module_QTL_signals$SNP_ID) )

# how many target genes for uni SNPs
table(battle_QTL_snp$n_gene)

# how many target genes for uni SNPs that are also multi SNPs
table(battle_QTL_snp %>% filter(if_rep) %>% select(n_gene))

# analyze uni SNPs that are not identified as multi SNPs
battle_QTL_snp %>%
  filter(!if_rep) %>%
  separate(snp, c('chr', 'pos'), ":", remove = FALSE) %>%
  left_join(battle_QTL %>% select(c(GENE_NAME, snp, LOG_PVAL)), by = "snp" ) %>%
  left_join(coexp_module, by = c("GENE_NAME" = "gene") ) %>%
  arrange(desc(n_gene)) %>%
  View()


### Plot the uni SNPs and those are also multi SNPs
df_bar <- battle_QTL_snp %>%
  group_by(n_gene) %>%
  summarise(n_SNP_uni = n(), n_SNP_multi = sum(if_rep) ) %>%
  pivot_longer(c("n_SNP_uni", "n_SNP_multi"), "SNP_type", values_to = "n_SNP")
df_bar$SNP_type <- factor(df_bar$SNP_type,
                          levels = c("n_SNP_uni", "n_SNP_multi"),
                          labels = c("Univariate", "Trans-PCO") )


fig_bar_prop <- ggplot(df_bar, aes(x = factor(n_gene), fill = SNP_type)) +
  geom_bar(aes(y = n_SNP), stat = "identity", position = "dodge2") +
  labs(x = "Number of Trans- Target Genes", y = "Number of Trans- SNPs", fill = "Method") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Paired")[c(1:2)]) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(color = "black", linetype = "dashed"),
        legend.key.size= unit(0.5, "cm"),
        axis.line = element_line(colour="black", size = 0.7),
        plot.margin=unit(c(10,5,5,5),"mm"),
        axis.text.x = element_text(colour = "black", angle = 0, hjust=1, vjust = 1, size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.y = element_text(angle=90,vjust =2, size=14),
        axis.title.x = element_text(vjust = -0.2, size=14),
        axis.title.y.right = element_text(angle = 90) )


### save figure and plotting object for further editing
saveRDS(fig_bar_prop, 'fig_compare_uni.rds')
ggsave("fig_compare_uni.png", fig_bar_prop,
       width = 5, height = 4)
