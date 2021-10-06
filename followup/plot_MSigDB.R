rm(list = ls())
require(ggplot2)
require(tidyr)


## Input files
file_gene_set = "h.all.v7.4.symbols.gmt"

## read file
gene_set = qusage::read.gmt(file_gene_set)

## gene set size
gene_set_N = length(gene_set)
gene_set_size = sort(sapply(gene_set, function(x) length(x)))
df_gene_set = data.frame("name" = names(gene_set_size),
                         "size" = gene_set_size,
                         stringsAsFactors = FALSE, check.names = FALSE)
df_gene_set$name = factor(df_gene_set$name, levels = df_gene_set$name)

## plot size
fig_bar <- ggplot(df_gene_set, aes(x = size, y = name)) +
  geom_bar(aes(fill = size), stat="identity") +
  labs(x = "Size", y = "Gene Set", fill = NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43")
ggsave(paste0("geneSet_size_", file_gene_set, ".png"), fig_bar)


## Overlapped genes among gene sets
overlap_mat = matrix(0, ncol = gene_set_N, nrow = gene_set_N,
                     dimnames = list(df_gene_set$name, df_gene_set$name))
for(i in 1:(gene_set_N-1)){
  for(j in (i+1):gene_set_N){
    overlap_mat[i, j] = sum(gene_set[[df_gene_set$name[i]]] %in% gene_set[[df_gene_set$name[j]]])
  }
}


## data frame for plotting
df_cast = cbind("geneSet1" = rownames(overlap_mat), as.data.frame(overlap_mat, check.names = FALSE, stringsAsFactors = FALSE), stringsAsFactors = FALSE) %>%
  pivot_longer(!geneSet1, names_to = "geneSet2", values_to = "overlap")
df_cast$geneSet1 = factor(df_cast$geneSet1, levels =  df_gene_set$name)
df_cast$geneSet2 = factor(df_cast$geneSet2, levels = df_gene_set$name )

### tile plot
fig_tile <- ggplot(df_cast, aes(geneSet2, geneSet1)) +
  geom_tile(aes(fill = overlap)) +
  scale_fill_gradient(low="white", high="grey2") +
  geom_text(data = df_cast[df_cast$overlap>10, ], aes(label = overlap), color = "tomato3", size = 2) +
  labs(x = "Gene set", y = "Gene set", fill = "#overlap genes") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), axis.text.y = element_text(size = 5))
ggsave(paste0("geneSet_overlap_", file_gene_set, ".png"), fig_tile, width = 10, height = 7)
system(paste0("bash ~/imgcat ", "geneSet_overlap_", file_gene_set, ".png"))

