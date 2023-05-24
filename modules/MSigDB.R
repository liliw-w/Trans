rm(list = ls())
require(ggplot2)

## Input files
file_gene_set = "h.all.v7.4.symbols.gmt"
file.ex.var.regressed = "~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds"
file.ex = "/project2/xuanyao/llw/DGN_PCO.lambda.01_real/data/ex.rds"
file_large_ex.var.regressed = "~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds"
out.coexp.module = "result/coexp.module.rds"
out_ex.var.regressed = "/scratch/midway2/liliw1/MODULES/MSigDB/ex.var.regressed.rds"

## read file
gene_set = qusage::read.gmt(file_gene_set)
datExpr = readRDS(file.ex.var.regressed)
ex = readRDS(file.ex)

## genes in datasets
gene_set_N = length(gene_set)
genes_msigDB = unique(unlist(gene_set))
genes_data = colnames(datExpr)
genes_origdata = sapply(colnames(ex), function(x) strsplit(x, "\\|")[[1]][1])

## gene sets with genes included in data
gene_set_in_data = lapply(gene_set, function(x) x[x %in% genes_data])
gene_set_size = sort(sapply(gene_set_in_data, function(x) length(x)), decreasing = TRUE)
gene_set_in_data = gene_set_in_data[names(gene_set_size)]

## convert gene sets to gene modules format, as input to Trans-PCO pipleine
out = list("moduleLabels" = setNames(rep(1:gene_set_N, gene_set_size), unlist(gene_set_in_data)))


## save and print
saveRDS(out, out.coexp.module)
cat("Out of", length(genes_msigDB), "MSigDB genes,",
    sum(genes_msigDB %in% genes_data), "are included in the dataset with", length(genes_data), "genes,",
    "while", sum(genes_msigDB %in% genes_origdata), "are included in the original data with", length(genes_origdata), "genes.",
    "\n")

## extract the ex.var.regressed with MSigDB genes from previous files
large_ex.var.regressed = readRDS(file_large_ex.var.regressed)
ex.var.regressed = large_ex.var.regressed[, names(out$moduleLabels)]
saveRDS(ex.var.regressed, out_ex.var.regressed)
