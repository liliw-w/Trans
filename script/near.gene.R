rm(list = ls())
parameter = commandArgs(trailingOnly = T)
file.qtl = parameter[1]
file.gene.meta = parameter[2]
file.gene = parameter[3]

file.qtl = 'postanalysis/indep.signals.chr.module.txt'
file.gene.meta = '/project2/xuanyao/llw/DGN/data/gene.meta.txt'
file.gene = 'postanalysis/near.gene.chr.module.indep.txt'

qtl = read.table(file.qtl,
                 sep = "\t", header = F, row.names = NULL,
                 stringsAsFactors = F)
qtl.df = data.frame("chr" = paste0("chr", sapply(qtl[, 1], function(x) strsplit(x, ":")[[1]][1])),
                 "pos" = as.numeric(sapply(qtl[, 1], function(x) strsplit(x, ":")[[1]][2])),
                 stringsAsFactors = F)
rownames(qtl.df) = qtl[, 1]

gene.meta = read.table(file.gene.meta,
                       sep ="\t", header=TRUE, row.names = NULL,
                       stringsAsFactors=F,  check.names = FALSE)
rownames(gene.meta) = gene.meta[, "gene"]


gene = lapply(1:nrow(qtl.df), function(x){
  snp = qtl.df[x, ]
  gene.near = which(gene.meta$chr==snp$chr & abs(gene.meta$start-snp$pos)<5*10^5)
  gene.nearest = gene.near[which.min(abs(gene.meta$start[gene.near]-snp$pos))]
  gene.meta[gene.nearest, "gene"]
  }
)

out = sapply(unique(unlist(gene)), function(x) strsplit(x, "|", fixed = T)[[1]][1])
write.table(out, file.gene, quote = F, sep = "\t", row.names = F, col.names = F)
