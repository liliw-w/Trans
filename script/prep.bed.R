#rm(list = ls())
#parameter = commandArgs(trailingOnly = T)
#dir.data = parameter[1]
#data.type = parameter[2]

data.type = snakemake@params[['data_type']]
input1 = snakemake@input[['file_sample_name']]
input2 = snakemake@input[['file_ex']]
input3 = snakemake@input[['file_covariates']]
input4 = snakemake@input[['file_gene_meta']]
input5 = snakemake@input[['file_coexp_module']]
file.covariates.null = snakemake@output[['file_covariates_null']]
file.expression = snakemake@output[['file_expression']]

str(file.expression)
print(file.expression)
print(file.expression[[1]])

options(stringsAsFactors = FALSE)

if(!data.type %in% c("obs", "null")){stop("Please specify a valid data.type: obs, null.\n")}

# sample name
sample.name = read.table(input1,
                         row.names = NULL, header=F,
                         stringsAsFactors=F, check.names = FALSE)

# expression
load(input2)
if(data.type == "null"){
  perm = sample(nrow(ex))
  ex = ex[perm, ]
  
  cov_all = read.table(input3,
                       sep = "\t", header = TRUE, row.names = 1,
                       stringsAsFactors = FALSE, check.names = FALSE)
  cov_all.perm = cov_all[, perm]; colnames(cov_all.perm) = colnames(cov_all)
  write.table(cov_all.perm, file.covariates.null, sep = "\t", quote = FALSE)
}
colnames(ex) = gnames; rownames(ex) = sample.name[, 1]

# gene position
gene.meta = as.matrix(read.table(input4,
                                 sep ="\t", header=TRUE, row.names = NULL,
                                 stringsAsFactors=F,  check.names = FALSE))
rownames(gene.meta) = gene.meta[, "gene"]
gene.meta = gene.meta[, c("chr", "start", "end", "gene")]
gene.meta[, "chr"] = paste0(gene.meta[, "chr"], "NA")
colnames(gene.meta)[1] = "#chr"

# module
module.info = readRDS(input5)
Nmodule = max(module.info$moduleLabels)
for(k in 1:Nmodule){
  gene.in.module = names(module.info$moduleLabels)[module.info$moduleLabels==k]
  res = cbind(gene.meta[gene.in.module, ], t(ex[, gene.in.module]))
  
  write.table(res, gzfile(file.expression[[k]]),
              sep = "\t", row.name = F, quote = F)
  
  nGene = nrow(res)
  nBatch = nGene %/% 100; nLeft = nGene %% 100
  if(nBatch > 0){
    for(i in 1:nBatch){
      write.table(res[(i*100-99):(i*100), ],
                  gzfile(paste0(file.expression[[k]], ".", i,".bed.gz")),
                  sep = "\t", row.name = F, col.names = T, quote = F)
    }
    if(nLeft != 0){
      write.table(res[(nBatch*100+1):nrow(res), ],
                  gzfile(paste0(file.expression[[k]], ".", (nBatch+1),".bed.gz")),
                  sep = "\t", row.name = F, col.names = T, quote = F)
    }
  }
}
