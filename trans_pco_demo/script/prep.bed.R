##############################################
########### Prepare bed file of a module to run tensorqtl ###########
##############################################
# load packages -----
library(data.table)
options(stringsAsFactors = FALSE)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    150
    data/ex.rds
    result/gene.meta.txt
    result/coexp.module.rds
    result/expression.module150.bed.gz
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

module <- as.numeric(args[1])
file.ex <- args[2]
file.gene.meta <- args[3]
file.coexp.module <- args[4]


## output -----
file.expression <- args[5]


# read files -----
ex <- readRDS(file.ex)
gene.meta <- fread(file.gene.meta, sep ="\t", header=TRUE)
coexp.module <- readRDS(file.coexp.module)


# gene position -----
gene.meta <- gene.meta[!duplicated(gene.meta$gene), ]
gene.meta <- gene.meta[, c("chr", "start", "end", "gene")]
gene.meta$chr <- paste0(gene.meta$chr, "NA")
colnames(gene.meta)[1] <- "#chr"
rownames(gene.meta) <- gene.meta$gene

Nmodule <- max(coexp.module$moduleLabels)


# format bed file for the module -----
gene.in.module <- names(coexp.module$moduleLabels)[coexp.module$moduleLabels == module]
res <- cbind(gene.meta[match(gene.in.module, gene.meta$gene), ], t(ex[, gene.in.module]))

fwrite(res, file.expression, sep = "\t", row.names = FALSE, col.names = TRUE)


## separate large modules in batches -----
nGene <- nrow(res)
nBatch <- nGene %/% 100
nLeft <- nGene %% 100
if(nBatch > 0){
  for(i in 1:nBatch){
    fwrite(
      res[(i*100-99):(i*100), ],
      paste0(file.expression, ".", i,".bed.gz"),
      sep = "\t", row.names = FALSE, col.names = TRUE
    )
  }
  if(nLeft != 0){
    fwrite(
      res[(nBatch*100+1):nrow(res), ],
      paste0(file.expression, ".", (nBatch+1),".bed.gz"),
      sep = "\t", row.names = FALSE, col.names = TRUE
    )
  }
}

