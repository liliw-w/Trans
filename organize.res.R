rm(list = ls())
require(data.table)

parameter = commandArgs(trailingOnly = T)
if(length(parameter) == 0){
  parameter = c('gencode.v19.annotation.table.TSS.txt',
                'FDR/signals.chr.module.perm10.txt',
                'postanalysis/LD.prun.in.chr.module.perm10.txt',
                "signals.nearset.gene.txt",
                "signals.nearset.gene.func.txt")
}

file_meta_TSS = parameter[1]
file_signals = parameter[2]
file_uniq_signal = parameter[3]
outfile_nearest_gene = parameter[4]
outfile_nearest_funcgene = parameter[5]

meta=as.data.frame(fread(file_meta_TSS, header=T))

snp=read.table(file_uniq_signal, header=F, stringsAsFactors=F, col.names='SNP')
snp=cbind(snp, "chr"=sapply(snp$SNP, function(x) paste0("chr", strsplit(x, ":")[[1]][1])), "pos"=sapply(snp$SNP, function(x) as.integer(strsplit(x, ":")[[1]][2])), stringsAsFactors=FALSE)

signals = as.data.frame(fread(file_signals, header=FALSE, col.names=c("signal", "p", "q")))
signals = cbind(signals,
                "SNP" = sapply(signals$signal, function(x) paste0(strsplit(x, ":")[[1]][-1], collapse=":")),
                "module" = sapply(signals$signal, function(x) strsplit(x, ":")[[1]][1]),
                stringsAsFactors=FALSE)


snp_gene = do.call(rbind, apply(snp, 1, function(x){
  ind = meta$Chromosome == as.character(x['chr'])

  intron = sign(meta[ind, "Start"] - as.integer(x['pos'])) * sign(meta[ind, "End"] - as.integer(x['pos']))
  if(any(intron == -1)){
    ind_mind = which(intron == -1)
    mind = abs(meta[ind, ][ind_mind, "TSS"] - as.integer(x['pos']))
    location = "within"
  }else{
    ind_mind = which.min(abs(meta[ind, "TSS"] - as.integer(x['pos'])))
    mind = min(abs(meta[ind, "TSS"] - as.integer(x['pos'])))
    location = "out"
  }

  ind_sig = signals[, "SNP"] %in% x["SNP"]

  if(any(mind<(5e+5))){
    return(do.call(rbind, lapply(which(mind<5e+5), function(y) cbind("SNP" = x["SNP"], "location" = location, "distance_TSS" = as.integer(mind)[y],
                                                              "minp" = min(signals[ind_sig, "p"]), "minq" = min(signals[ind_sig, "q"]),
                                                              "module"= paste0(signals[ind_sig, "module"], collapse=";"),
                                                              meta[ind, ][ind_mind[y], ], stringsAsFactors=FALSE, row.names = NULL)) ))
  }else{
    res = cbind("SNP" = x["SNP"], "location" = location, "distance_TSS" = as.integer(mind),
                "minp" = min(signals[ind_sig, "p"]), "minq" = min(signals[ind_sig, "q"]),
                "module"= paste0(signals[ind_sig, "module"], collapse=";"),
                meta[ind, ][ind_mind, ], stringsAsFactors=FALSE, row.names = NULL)
    res[c(-1, -4, -5, -6)] = NA
    return(res)
  }
}
))
cat(sum(is.na(snp_gene$Chromosome)), "SNPs have no nearset gene.", "\n")
cat(length(unique(snp_gene$GeneSymbol)), "nearest genes in total.", "\n")

fwrite(snp_gene[order(snp_gene$minp), ], outfile_nearest_gene, quote=FALSE, sep='\t', na = NA)
cat("File", outfile_nearest_gene, "completed!", "\n")


meta_func = meta[meta$Class %in% c('protein_coding', 'lincRNA'), ]
snp_gene_func = do.call(rbind, apply(snp, 1, function(x){
  ind = meta_func$Chromosome == as.character(x['chr'])

  intron = sign(meta_func[ind, "Start"] - as.integer(x['pos'])) * sign(meta_func[ind, "End"] - as.integer(x['pos']))
  if(any(intron == -1)){
    ind_mind = which(intron == -1)
    mind = abs(meta_func[ind, ][ind_mind, "TSS"] - as.integer(x['pos']))
    location = "within"
  }else{
    ind_mind = which.min(abs(meta_func[ind, "TSS"] - as.integer(x['pos'])))
    mind = min(abs(meta_func[ind, "TSS"] - as.integer(x['pos'])))
    location = "out"
  }

  ind_sig = signals[, "SNP"] %in% x["SNP"]

  if(any(mind<(5e+5))){
    return(do.call(rbind, lapply(which(mind<5e+5), function(y) cbind("SNP" = x["SNP"], "location" = location, "distance_TSS" = as.integer(mind)[y],
                                                                     "minp" = min(signals[ind_sig, "p"]), "minq" = min(signals[ind_sig, "q"]),
                                                                     "module"= paste0(signals[ind_sig, "module"], collapse=";"),
                                                                     meta_func[ind, ][ind_mind[y], ], stringsAsFactors=FALSE, row.names = NULL)) ))
  }else{
    res = cbind("SNP" = x["SNP"], "location" = location, "distance_TSS" = as.integer(mind),
                "minp" = min(signals[ind_sig, "p"]), "minq" = min(signals[ind_sig, "q"]),
                "module"= paste0(signals[ind_sig, "module"], collapse=";"),
                meta[ind, ][ind_mind, ], stringsAsFactors=FALSE, row.names = NULL)
    res[c(-1, -4, -5, -6)] = NA
    return(res)
  }
}
))
cat(sum(is.na(snp_gene_func$Chromosome)), "SNPs have no nearset gene.", "\n")
cat(length(unique(snp_gene_func$GeneSymbol)), "nearest genes in total.", "\n")

fwrite(snp_gene_func[order(snp_gene_func$minp), ], outfile_nearest_funcgene, quote=FALSE, sep='\t', na = NA)
cat("File", outfile_nearest_funcgene, "completed!", "\n")

