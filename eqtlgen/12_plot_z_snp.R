##############################################
########### plot z-scores of (one SNP, modules) in DGN and eQTLGen ###########
########### for SNPs are eQTLGen specific signals not in DGN ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)

source('~/Trans/plot/theme_my_pub.R')


# paras and I/O -----
module <- 54
chr <- 3
snp_list <- c("3:47262246", "3:49731861", "3:56849749")

file_z_eqtlgen <- paste0("/project2/xuanyao/llw/eQTLGen/z_dgn_166_module/z.module", module, ".chr", chr, ".txt.gz")
file_z_dgn <- paste0("/project2/xuanyao/llw/DGN_no_filter_on_mappability/z/z.module", module, ".chr", chr, ".txt.gz")

file_snp_meta <- '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'
file_gene_meta <- '/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt'


# read files -----
z_eqtlgen <- fread(file_z_eqtlgen)
z_dgn <- fread(file_z_dgn)

snp_meta <- fread(file_snp_meta)
gene_meta <- fread(file_gene_meta)


# convert snp ids -----
z_eqtlgen <- z_eqtlgen %>%
  left_join(snp_meta %>% select(SNP, meta), by = c("snp" = "SNP")) %>%
  relocate(meta, .after = "snp") %>%
  select(-snp) %>%
  rename("snp" = meta)

# convert gene ids -----
gene_meta <- gene_meta %>% separate(Geneid, into = c("gene", NA), sep = "[.]")
colnames(z_eqtlgen)[-1] <- gene_meta[match(colnames(z_eqtlgen)[-1], gene_meta$gene), GeneSymbol]

z_dgn <- z_dgn %>% select(colnames(z_eqtlgen))


# select snps of interest to draw -----
plt_df <- rbind(
  z_dgn %>%
    filter(snp %in% snp_list) %>%
    pivot_longer(-snp, names_to = "gene", values_to = "z") %>%
    mutate("type" = "DGN"),
  
  z_eqtlgen %>%
    filter(snp %in% snp_list) %>%
    pivot_longer(-snp, names_to = "gene", values_to = "z") %>%
    mutate("type" = "eQTLGen")
)


# Violin plot of SNPs v.s. their z-scores in two datasets -----
base_plt <- ggplot(plt_df, aes(x = snp, y = z)) +
  geom_violin(aes(fill = type)) +
  geom_point(
    aes(color = type),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 1),
    shape = 16, size = 1, alpha = 1
  ) +
  labs(x = NULL, y = "Z-score", color = "Dataset", title = paste("Module", module))


base_plt +
  scale_colour_manual(
    breaks = c("DGN", "eQTLGen"),
    values = c("eQTLGen" = "#0028a1", "DGN" = "#85192d"),
    guide = guide_legend(override.aes = list(size = 2, alpha = 1))
  ) +
  scale_fill_manual(
    breaks = c("DGN", "eQTLGen"),
    values = c("eQTLGen" = "#d6dcef", "DGN" = "#e1c7cc"),
    guide = "none"
  ) +
  theme_my_pub(
    legend.position = "right",
    axis.text.x.angle = 60,
    axis.title.size = 14
  ) +
  theme(
    panel.grid.major.y = element_line(linetype = "dashed", color = "#999999"),
    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text.x = element_text(size = 10, hjust = 1, vjust = 1),
    
    plot.title = element_text(vjust = -1)
  )

ggsave(paste0("plot/z_eqtlgen_dgn_M", module, "_chr", chr, ".pdf"), height = 5, width = 6)

