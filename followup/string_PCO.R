##############################################################################
######### Prepare coexp.module and ex.var.regressed files to run PCO #########
##############################################################################
rm(list = ls())

## Input files
file_gene_set = "/scratch/midway2/liliw1/MODULES/STRING/module.rds"
file.ex = "/project2/xuanyao/llw/DGN_PCO.lambda.01_real/data/ex.rds"
file_large_ex.var.regressed = "~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds"
out.coexp.module = "result/coexp.module.rds"
out_ex.var.regressed = "result/ex.var.regressed.rds"

## read file
gene_set = readRDS(file_gene_set)
large_ex.var.regressed = readRDS(file_large_ex.var.regressed)
ex = readRDS(file.ex)

## genes in datasets
genes_data = colnames(large_ex.var.regressed)
genes_origdata = sapply(colnames(ex), function(x) strsplit(x, "\\|")[[1]][1])

## convert gene sets to gene modules format, as input to Trans-PCO pipleine
out = list("moduleLabels" = sort(gene_set[names(gene_set) %in% genes_data]) )

## save and print
saveRDS(out, out.coexp.module)
cat("Out of", length(gene_set), "MSigDB genes,",
    sum(names(gene_set) %in% genes_data), "are included in the dataset with", length(genes_data), "genes,",
    "while", sum(names(gene_set) %in% genes_origdata), "are included in the original data with", length(genes_origdata), "genes.",
    "\n")

## extract the ex.var.regressed with MSigDB genes from previous files
ex.var.regressed = large_ex.var.regressed[, names(out$moduleLabels)]
saveRDS(ex.var.regressed, out_ex.var.regressed)
