dir_data = '/home/liliw1/xuanyao_data/DGN/genotype/'

load(paste0(dir_data, 'ex.rdata'))
gene_meta = read.table(paste0(dir_data, 'gene.meta.txt'),
               sep = '\t', row.names = NULL, header = T,
               stringsAsFactors = F)



merged = merge(data.frame("gene" = gnames, stringsAsFactors = F), gene_meta, by.x=1, by.y=1)

if(sum(!(gnames %in% merged$gene))){
  ex = ex[, gnames %in% merged$gene]
  gnames = gnames[gnames %in% merged$gene]
}

save(ex, gnames, file = paste0(dir_data, 'ex.rdata'))
