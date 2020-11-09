# dir.data is the directory where gene meta files are stored.
# Files should include:
#  1. gene.meta.txt: gene meta file, with header: gene, chr, start, end
#  2. pseudogenes.txt: two columns, separated by \t.
#  3. mappability.txt: two columns, separated by \t.

rm(list = ls())

#parameter = commandArgs(trailingOnly = T)
#dir.data = parameter[1]

options(stringsAsFactors = FALSE)

### Data preparation
datExpr = readRDS(snakemake@input[['file_ex_var_regressed']])
gene.name = sapply(colnames(datExpr), function(x) strsplit(x, "\\|")[[1]][1])


# remove genes on chromosome X & Y
gene.meta = read.table(snakemake@input[['file_gene_meta']],
                       header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rownames(gene.meta) = gene.meta$gene
ind_chrX = gene.meta[colnames(datExpr), "chr"] == "chrX" | gene.meta[colnames(datExpr), "chr"] == "chrY"

# remove pseudogenes
pseudogenes = read.table(snakemake@input[['file_pseudogenes']], sep ="\t", stringsAsFactors = F)
ind_pseudo = gene.name %in% unique(pseudogenes[, 1])

# remove low mappability genes
mappability = read.table(snakemake@input[['file_mappability']], sep ="\t", stringsAsFactors = F)
ind_map = gene.name %in% unique(mappability[, 2])

# remove cross-mappable genes
cross_mappability = readRDS(snakemake@input[['file_cross_mappable_genes']])
cross_map_genes = unique(names(cross_mappability[!is.na(cross_mappability)]))
ind_cross_map = gene.name %in% cross_map_genes

ind_remove = ind_chrX | ind_pseudo | ind_map | ind_cross_map

sum(ind_chrX)
sum(ind_pseudo)
sum(ind_map)
sum(ind_cross_map)
sum(ind_remove)


write.table(as.matrix(data.frame("gene" = colnames(datExpr),
                                 "gene.name" = gene.name,
                                 "ind_remove" = ind_remove,
                                 "ind_chrX" = ind_chrX,
                                 "pseudo" = ind_pseudo,
                                 "lowmap" = ind_map,
                                 "crossmap" = ind_cross_map)),
            snakemake@output[["file_genes_rm_info"]],
            sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
