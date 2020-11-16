rm(list = ls())

parameter = commandArgs(trailingOnly = T)
signals_file = parameter[1]
sig_uniq = parameter[2]


signals = read.table(signals_file,
                     sep = "\t", header = FALSE, row.names = 1,
                     stringsAsFactors = FALSE)

module.snp = rownames(signals)
snp = sort(unique(sapply(strsplit(module.snp, ":"), function(x) paste(x[-1], collapse = ":"))))

write.table(snp, file = sig_uniq, quote = FALSE, row.names = FALSE, col.names = FALSE)
