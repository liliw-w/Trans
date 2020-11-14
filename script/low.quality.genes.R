#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#dir.data = parameter[1]

input1 = snakemake@input[['file_ex_var_regressed']]
input2 = snakemake@input[['file_gene_meta']]
input3 = snakemake@input[['file_pseudogenes']]
input4 = snakemake@input[['file_mappability']]
input5 = snakemake@input[['file_cross_mappable_genes']]
output1 = snakemake@output[["file_genes_rm_info"]]

options(stringsAsFactors = FALSE)

### Data preparation
datExpr = readRDS(input1)
gene.name = sapply(colnames(datExpr), function(x) strsplit(x, "\\|")[[1]][1])


# remove genes on chromosome X & Y
gene.meta = read.table(input2, sep ="\t",
                       header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rownames(gene.meta) = gene.meta$gene

chr.auto = paste0("chr", 1:22)
ind_chrXY = !(gene.meta[colnames(datExpr), "chr"] %in% chr.auto)

# remove pseudogenes
pseudogenes = read.table(input3, sep ="\t", stringsAsFactors = F)
ind_pseudo = gene.name %in% unique(pseudogenes[, 1])

# remove low mappability genes
mappability = read.table(input4, sep ="\t", stringsAsFactors = F)
ind_map = gene.name %in% unique(mappability[, 2])

# remove cross-mappable genes
cross_mappability = readRDS(input5)
cross_map_genes = unique(names(cross_mappability[!is.na(cross_mappability)]))
ind_cross_map = gene.name %in% cross_map_genes

ind_remove = ind_chrXY | ind_pseudo | ind_map | ind_cross_map

cat("All genes:", ncol(datExpr), "\n")
cat("XY genes:", sum(ind_chrXY), "\n")
cat("pseodo genes:", sum(ind_pseudo), "\n")
cat("low map genes:", sum(ind_map), "\n")
cat("cross map genes:", sum(ind_cross_map), "\n")
cat("Remove genes:", sum(ind_remove), "\n")
cat("Left genes:", ncol(datExpr)-sum(ind_remove), "\n")


write.table(as.matrix(data.frame("gene" = colnames(datExpr),
                                 "gene.name" = gene.name,
                                 "ind_remove" = ind_remove,
                                 "ind_chrXY" = ind_chrXY,
                                 "pseudo" = ind_pseudo,
                                 "lowmap" = ind_map,
                                 "crossmap" = ind_cross_map)),
            output1,
            sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
