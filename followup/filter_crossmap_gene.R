rm(list = ls())
library(data.table)

## files
file_crossGene_gene = '/scratch/midway2/liliw1/DGN_cross_map_filter/crossGene_gene.rds'
file_cross_mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"
file_out = "crossmap_gene_2000.txt"
thre_num_crossmap = 2000

## read data
crossGene_gene = readRDS(file_crossGene_gene)
cross_map = fread(file_cross_mappability)
thre_crossmap_score = quantile(cross_map$V3, 0.9)

## filter criteria 1: genes with too many cross-map genes
rem_gene1 = names(crossGene_gene[crossGene_gene > thre_num_crossmap]) #4314

## filter criteria 2: gene pairs with too high cross-map score
#rem_gene2 = unique(cross_map[cross_map$V3 > thre_crossmap_score, c(V1, V2)])
#rem_gene2 = rem_gene2[rem_gene2 %in% names(crossGene_gene)]


## write out
#str(rem_gene1)
#str(rem_gene2)
#str(unique(c(rem_gene1, rem_gene2)))
rem_gene = rem_gene1
fwrite(as.data.frame(rem_gene1), file = file_out, quote = FALSE, sep = "\t")

#quantile(cross_map$V3, c(0.75, 0.8, 0.9) )
#75%  80%  90%
#17.0 20.5 33.0
