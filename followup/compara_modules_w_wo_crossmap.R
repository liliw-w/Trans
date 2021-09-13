rm(list = ls())
require(tidyr)
require(ggplot2)

## File I/O
file.coexp.module.noCross = "/scratch/midway2/liliw1/DGN_PCO.lambda.01_real/result/coexp.module.rds"
file.coexp.module.Cross = "/scratch/midway2/liliw1/DGN_cross_map_filter/result_1000/coexp.module.rds"
fig_file = "plot/overlap_gene_across_modules.png"

## read data
coexp_module_noCross = readRDS(file.coexp.module.noCross)$moduleLabels
coexp_module_Cross = readRDS(file.coexp.module.Cross)$moduleLabels

## parameters and result storage
Nmodule_Cross = max(coexp_module_Cross); Nmodule_noCross = max(coexp_module_noCross)
mat_overlap = matrix(nrow = Nmodule_Cross, ncol = Nmodule_noCross,
                     dimnames = list(paste0("cM", 1:Nmodule_Cross), paste0("ncM", 1:Nmodule_noCross)) )
mat_overlap = cbind(mat_overlap, "size" = as.numeric(table(coexp_module_Cross)[as.character(1:Nmodule_Cross)]) )
mat_overlap = rbind(mat_overlap, "size" = c(as.numeric(table(coexp_module_noCross)[as.character(1:Nmodule_noCross)]), 0) )

## calculate overlap genes among modules from two different settings
for(i in 1:Nmodule_Cross){
  gene_within_module = names(coexp_module_Cross[coexp_module_Cross == i])
  mat_overlap[i, -ncol(mat_overlap)] = sapply(1:Nmodule_noCross, function(x){
    gene_noCross = names(coexp_module_noCross[coexp_module_noCross == x])
    sum(gene_within_module %in% gene_noCross)
  } )
}

## only leave the top modules
ind_top_module_Cross = c(nrow(mat_overlap), 1:max(as.numeric(names(which(table(coexp_module_Cross)>100)))) )
most_overlap_noCross = c("size", unique(colnames(mat_overlap)[apply(mat_overlap[, -ncol(mat_overlap)], 1, which.max)][ind_top_module_Cross]))
mat_overlap_top = mat_overlap[ind_top_module_Cross, most_overlap_noCross]

## Tile plot: overlap genes across two module networks
### prepare data
df_cast = cbind("Cross" = rownames(mat_overlap_top), as.data.frame(mat_overlap_top, check.names = FALSE, stringsAsFactors = FALSE), stringsAsFactors = FALSE) %>%
  pivot_longer(!Cross, names_to = "noCross", values_to = "overlap")
df_cast$Cross = factor(df_cast$Cross, levels = c("size", rownames(mat_overlap[-nrow(mat_overlap),])) )
df_cast$noCross = factor(df_cast$noCross, levels = c("size", colnames(mat_overlap[,-ncol(mat_overlap)])) )

### tile plot
fig_tile <- ggplot(df_cast, aes(noCross, Cross)) +
  geom_tile(aes(fill = overlap)) +
  scale_fill_gradient(low="white", high="grey2") +
  geom_text(data = df_cast[df_cast$overlap>10, ], aes(label = overlap), color = "tomato3") +
  labs(x = "Modules without cross-map genes", y = "Modules with cross-map genes", fill = "#overlap genes") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top")

ggsave(fig_file, fig_tile)
