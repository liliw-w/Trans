##############################################
########### Calculate pvalue of association between a snp and gene module by PCO ###########
##############################################
# load packages -----
library(data.table)


# I/O & paras -----
if(interactive()){
  args <- scan(
    text = '
    result/ex.var.regressed.rds
    result/gene.meta.txt
    result/coexp.module.rds
    z/z.module150.chr20.txt.gz
    20
    150
    p/p.module150.chr20.rds
    ',
    what = 'character'
  )
} else{
  args <- commandArgs(trailingOnly = TRUE)
}

file.ex.var.regressed <- args[1]
file.gene.meta <- args[2]
file.coexp.module <- args[3]
file.z <- args[4]
chr <- as.numeric(args[5])
module <- as.numeric(args[6])


## output -----
file.p <- args[7]



# read files -----
datExpr <- readRDS(file.ex.var.regressed)
gene.meta <- read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module <- readRDS(file.coexp.module)
z.mat <- fread(file.z)
z.mat <- as.matrix(z.mat, rownames = TRUE)


# load PCO -----
source(paste0("script/ModifiedPCOMerged_acat.R"))
source(paste0("script/liu.R"))
source(paste0("script/liumod.R"))
source(paste0("script/davies.R"))
dyn.load(paste0("script/qfc.so"))


# select genes in trans for the target chromosome and only use these trans genes in the module -----
gene_in_cluster <- data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
gene_w_pos <- merge(gene_in_cluster, gene.meta)
gene_trans <- gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]


# run PCO -----
if(length(gene_trans) > 1){
  Sigma <- cor(datExpr[, gene_trans])
  z.mat_trans <- z.mat[, gene_trans]; rm(z.mat)
  
  p.all <- ModifiedPCOMerged_acat(Z.mat = z.mat_trans, Sigma = Sigma)
  
}else{
  p.all <- NULL
}

# save p -----
saveRDS(p.all, file = file.p)

