rm(list = ls())

parameter = commandArgs(trailingOnly = T)
sig = parameter[1]
sig.indp = parameter[2]
FDR.method = parameter[3]

library(ggplot2)
library(gridExtra)

signals = read.table(sig,
                       sep = "\t",
                       row.names = NULL, header = F, stringsAsFactors = F)
colnames(signals) = c("module.snp", "p", "q")
signals[signals$p == 0, "p"] = min(signals$p[signals$p!=0])/100
signals$snp = sapply(strsplit(signals[, 1], ":"), function(x) paste(x[-1], collapse = ":"))
signals$module = sapply(strsplit(signals[, 1], ":"), function(x) as.numeric(strsplit(x[1], "C")[[1]][2]))
signals$chr = sapply(strsplit(signals[, 1], ":"), function(x) as.numeric(strsplit(x[-1], ":")[[1]][1]))

indp.signals = read.table(sig.indp,
                      row.names = NULL, header = F, stringsAsFactors = F, col.names = "snp")
indp.signals$minp = sapply(indp.signals$snp, function(x) min(signals[signals$snp == x, "p"]))
indp.signals$chr = as.character(sapply(strsplit(indp.signals[, "snp"], ":"), function(x) x[1]))
indp.signals$chr = factor(indp.signals$chr, levels=1:22)

# histogram of #signals on chr
fig2 = ggplot(indp.signals, aes(chr)) + geom_histogram(stat="count")

# histogram of -logp
fig1 = ggplot(indp.signals, aes(x=-log10(minp))) + geom_histogram(binwidth=0.5) + coord_cartesian(xlim=c(0,20))

fig.all = list(fig1, fig2)
ggsave(paste0("plots/plot", FDR.method, ".png"), marrangeGrob(fig.all, nrow=1, ncol=2, top = NULL), width = 10, height = 5)
