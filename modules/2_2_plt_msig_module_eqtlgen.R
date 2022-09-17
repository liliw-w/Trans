##############################################
########### plot size of MSigDB modules with genes overlapped with eQTLGen ###########
##############################################
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggforce)
source('~/Trans/plot/theme_my_pub.R')



# I/O & paras -----
file_msigdb_module_df <- '/project2/xuanyao/llw/MODULES/MSigDB/h.all.v7.4.symbols.txt'
file_eqtlgen_genes <- "/project2/xuanyao/llw/eQTLGen/eqtlgen_genes.txt"
file_msigdb_dgn_module <- '/project2/xuanyao/llw/MODULES/MSigDB/result/coexp.module.rds'


## output -----
file_plt_module_size <- 'msigdb_module_eqtlgen.pdf'
file_plt_n_module <- 'nmodule_msigdb_gene.pdf'


# read files -----
msigdb_module <- fread(file_msigdb_module_df, header = TRUE)
eqtlgen_genes <- fread(file_eqtlgen_genes)
msigdb_dgn_module <- readRDS(file_msigdb_dgn_module)


# organize data -----
## overlap info with eqtlgen genes -----
msigdb_module <- msigdb_module %>%
  group_by(category) %>%
  mutate(module_size = n()) %>%
  ungroup() %>%
  mutate(if_eqtlgen = gene %in% !!eqtlgen_genes$GeneSymbol,
         anno = category,
         category = fct_infreq(category))



# plot module size and size of overlapped genes with eqtlgen -----
ggplot(msigdb_module, aes(category)) +
  geom_bar(aes(fill = "msigdb")) +
  geom_bar(data = filter(msigdb_module, if_eqtlgen), aes(fill = "eqtlgen")) +
  geom_bar(data = tmp, aes(fill = "dgn")) +
  labs(fill = "Data", x = "Module", y = "Size") +
  scale_fill_manual(
    breaks = c("msigdb", "eqtlgen", "dgn"),
    values = c("msigdb" = "#d9eaf3", "eqtlgen" = "#186191", "dgn" = "#8b1919"),
    labels = c("msigdb" = "MSigDB", "eqtlgen" = "Genes in eQTLGen", "dgn" = "Genes in DGN")
  ) +
  theme_my_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    panel.grid.major.y = element_line(linetype = "dashed"),
    legend.position = "bottom"
  )

ggsave(
  file_plt_module_size,
  width = 7, height = 5
)



# plot how many modules do every gene belongs to -----
ggplot(msigdb_module %>% group_by(gene) %>% summarise(n_module = n()),
       aes(n_module)) +
  geom_bar() +
  facet_zoom(xy = n_module > 3, zoom.size = 0.7) +
  labs(x = "Number of modules", y = "Gene count") +
  scale_x_continuous(breaks = seq(0, 100, 1)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file_plt_n_module,
  width = 6, height = 3
)

