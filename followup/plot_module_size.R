rm(list = ls())
require(ggplot2)
require(data.table)
require(dplyr)

# File name
file_coexp_module = "/scratch/midway2/liliw1/DGN_cross_map_filter/result/coexp.module.rds"
file_gene_annotation = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
file_cross_mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"

# Load data from files
coexp_module = readRDS(file_coexp_module)$moduleLabels
gene_annotation = fread(file_gene_annotation, header = TRUE)
cross_map = fread(file_cross_mappability)

# DF for each gene
df_module = data.frame("gene" = names(coexp_module),
                       "module" = coexp_module,
                       "if_cross_map" = names(coexp_module) %in% gene_annotation[match(unique(c(cross_map$V1, cross_map$V2)), gene_annotation$Geneid), GeneSymbol],
                       "num_cross_map_within_module" = NA,
                       check.names = FALSE, stringsAsFactors = FALSE)
### df_module %>% group_by(module) %>% summarise(module_size = n(), num_cross_map = sum(if_cross_map))

# Plot
fig <- ggplot(df_module) + geom_bar(aes(module, fill = if_cross_map)) +
  labs(x = "module", y = "size", fill = "if_cross_map") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top") +
  scale_fill_manual(values=c("dodgerblue2", "tomato3"))

ggsave("plot_module_size.png", fig)
system("bash ~/imgcat plot_module_size.png")


# Count genes with cross-mappable genes within each module
for (i in max(df_module$module):0) {
  K = sum(df_module$module==i)
  gene_names = df_module[df_module$module==i, "gene"]
  gene_names_anno = gene_annotation[match(gene_names, gene_annotation$GeneSymbol), Geneid]

  cross_within_module = matrix(nrow = K, ncol = K,
                               dimnames = list(gene_names_anno, gene_names_anno) )
  diag(cross_within_module) = FALSE
  for( j in 1:(K - 1) ){
    target = rownames(cross_within_module)[j]
    target_search = colnames(cross_within_module)[(j+1):K]

    target_cross_ind = cross_map$V1 %in% target | cross_map$V2 %in% target
    if(sum(target_cross_ind) == 0){
      cross_within_module[j, (j+1):K] = FALSE
      next
    }
    target_cross = cross_map[target_cross_ind, ]
    cross_within_module[j, (j+1):K] = target_search %in% c(target_cross$V1, target_cross$V2)
  }

  cross_within_module[lower.tri(cross_within_module)] = t(cross_within_module)[lower.tri(cross_within_module)]
  df_module[df_module$module==i, "num_cross_map_within_module"] = rowSums(cross_within_module)
  cat("Module: ", i, " is done!", "\n")
}
saveRDS(df_module, "df_module.rds")

fig_cross_within_module <- ggplot(df_module) + geom_bar(aes(module, fill = ifelse(num_cross_map_within_module==0, FALSE, TRUE))) +
  labs(x = "module", y = "size", fill = "if_cross_map_within_module") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top") +
  scale_fill_manual(values=c("dodgerblue2", "tomato3"))

ggsave("plot_module_size_cross_within_module.png", fig_cross_within_module)
system("bash ~/imgcat plot_module_size_cross_within_module.png")

