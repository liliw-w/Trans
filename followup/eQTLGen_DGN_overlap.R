# Nmodule: 75
rm(list = ls())
library(data.table)
library(ggplot2)
library(gridExtra)
qqplot.hist <- function(input, title){
  observed <- sort(input)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected+1))))
  df = data.frame(x = lexp, y = lobs, yy = observed)

  res = list()
  res[[1]] = ggplot(df, aes(x=x, y=y)) + geom_point(size = 0.2) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    theme(text=element_text(size=6)) +
    labs(title = title,
         x = "Expected (-logP)", y = "Observed (-logP)")
  res[[2]] = ggplot(df, aes(yy)) + geom_histogram(breaks = seq(0, 1, 0.01)) +
    theme(text=element_text(size=6)) +
    labs(title = title, x = "Observed (P)")

  return(res)
}

file_eQTLGen_snp = '/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt'
dir_p = "~/xuanyao_llw/eQTLGen_DGN_PCO.lambda.01/"
file.gene.meta = paste0(dir_p, "result/gene.meta.txt")
meta_used_snp = fread(file_eQTLGen_snp, header = TRUE)
file.p = list.files(path = paste0(dir_p, "p"), pattern = "^p.module.*chr.*rds", full.names = FALSE)

## Extract pvalus ##
p.obs = rbindlist(lapply(file.p, function(x)
{tmp_y=readRDS(paste0(dir_p, "p/", x));
res = tmp_y[names(tmp_y) %in% meta_used_snp$meta]
names(res) = meta_used_snp[match(names(res), meta_used_snp$meta), SNP]
saveRDS(res, paste0("p/", x));
print(x);
print(length(res));
as.data.table(setNames(res, paste0(strsplit(x, '.', fixed=T)[[1]][2], ":", names(res))), keep.rownames=T)
}))
setnames(p.obs, c('module_snp', 'p'))
p.obs$module = sapply(p.obs$module_snp, function(x) strsplit(x, ":")[[1]][1])
p.obs$snp_ID = sapply(p.obs$module_snp, function(x) strsplit(x, ":")[[1]][2])
p.obs$snp = meta_used_snp[match(p.obs$snp_ID, meta_used_snp$SNP), meta]
p.obs$snp_chr = sapply(p.obs$snp, function(x) strsplit(x, ":")[[1]][1])
p.obs = as.data.frame(p.obs[order(p.obs$p), ])
saveRDS(p.obs, "p.obs.all.rds")

## plot p value distribution and qq plots ##
fig.all = NULL
for (dir_pall in c("../eQTLGen_PCO.lambda.01/", "../eQTLGen_DGN_overlap/")) {
  p.obs = readRDS(paste0(dir_pall, "p.obs.all.rds"))

  fig.all = rbind(fig.all, qqplot.hist(p.obs$p, dir_p))
  cat(dir_p, "done!\n")
}
ggsave('qqplot.all.p.png',
       marrangeGrob(fig.all, nrow=length(fig.all)/2, ncol=2, top = NULL),
       height = length(fig.all)/2*(6/2), width = 6)

## Extract zscores ##
gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
file.z = list.files(path = paste0(dir_p, "z"), pattern = "^z.module.*chr.*txt.gz", full.names = FALSE)
for(x in file.z){
  z.eQTLGen = fread(paste0("~/scratch/eQTLGen_PCO.lambda.01/z/", x))
  z.mat = fread(paste0(dir_p, "z/", x))

  setnames(z.mat, c("snp", gene.meta[match(colnames(z.mat)[-1], gene.meta$gene), "GeneNameConv"]))
  ind_gene = colnames(z.mat) %in% colnames(z.eQTLGen)
  z.mat = z.mat[, ..ind_gene]
  z.mat = z.mat[z.mat$snp %in% meta_used_snp$meta, ]
  z.mat$snp = meta_used_snp[match(z.mat$snp, meta_used_snp$meta), SNP]

  order_gene = match(colnames(z.eQTLGen), colnames(z.mat))
  z.mat = z.mat[, ..order_gene]

  fwrite(z.mat, paste0("z/", x),
         sep = "\t", row.names = FALSE, col.names = TRUE)
  print(x)
}

module = 1; chr = 1
file_z = paste0("z.module", module, ".chr", chr, ".txt.gz")
file_p = paste0("p.module", module, ".chr", chr, ".rds")
z.eQTLGen = fread(paste0("~/scratch/eQTLGen_PCO.lambda.01/z/", file_z))
z.DGN.eQTLGen = fread(paste0("~/scratch/eQTLGen_DGN_overlap/z/", file_z))
p.eQTLGen = readRDS(paste0("~/scratch/eQTLGen_PCO.lambda.01/p/", file_p))
p.DGN.eQTLGen = readRDS(paste0("~/scratch/eQTLGen_DGN_overlap/p/", file_p))

z.eQTLGen = z.eQTLGen[z.eQTLGen$snp %in% z.DGN.eQTLGen$snp, ]
z.eQTLGen = cbind("p" = p.eQTLGen[match(z.eQTLGen$snp, names(p.eQTLGen))],
                  "type" = "eQTLGen",
                  z.eQTLGen)
z.eQTLGen = z.eQTLGen[order(z.eQTLGen$p), ]

z.DGN.eQTLGen = z.DGN.eQTLGen[match(z.eQTLGen$snp, z.DGN.eQTLGen$snp), ]
z.DGN.eQTLGen = cbind("p" = p.DGN.eQTLGen[match(z.DGN.eQTLGen$snp, names(p.DGN.eQTLGen))],
                      "type" = "DGN.eQTLGen",
                      z.DGN.eQTLGen)

z.eQTLGen = cbind("maxz" = apply(z.eQTLGen[, -(1:3)], 1, function(x) max(abs(x))),
                  z.eQTLGen)
z.DGN.eQTLGen = cbind("maxz" = apply(z.DGN.eQTLGen[, -(1:3)], 1, function(x) max(abs(x))),
                  z.DGN.eQTLGen)

png("tmp.png")
par(mfrow = c(1, 2))
plot(unlist(z.eQTLGen[1, -(1:4)]), col = "red", ylim = c(-4, 3))
points(unlist(z.DGN.eQTLGen[5, -(1:4)]), col = "black")
plot(unlist(z.eQTLGen[3, -(1:4)]), col = "red", pch = 18, ylim = c(-4, 3))
points(unlist(z.DGN.eQTLGen[3, -(1:4)]), col = "black", pch = 18)
dev.off()

"rs928391"
range(z.eQTLGen$p)
range(z.DGN.eQTLGen$p)


z.eQTLGen[1:3, 1:5]
z.DGN.eQTLGen[1:3, 1:5]

identical(z.DGN.eQTLGen$snp, z.eQTLGen$snp)
identical(colnames(z.DGN.eQTLGen), colnames(z.eQTLGen))


eQTLGen_signal = fread("2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz", header = TRUE)
str(unique(eQTLGen_signal$SNP))

require(data.table)

Comp = fread("/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz", header = TRUE)
Us_sig = fread("/project2/xuanyao/llw/eQTLGen_DGN_PCO.lambda.01/postanalysis/LD.prun.in.chr.module.perm10.txt", header = FALSE, col.names = "sig")
meta = fread("/project2/xuanyao/llw/eQTLGen_DGN_PCO.lambda.01/result/gene.meta.txt", header = TRUE)

Comp$snp = paste(Comp$SNPChr, Comp$SNPPos, sep=":")

str(unique(Comp$snp))
str(unique(Comp$GeneSymbol))

gene_us = meta[!meta$Remove, gene]
Comp = Comp[Comp$GeneSymbol %in% gene_us, ]

str(unique(Comp$snp))
str(unique(Comp$GeneSymbol))

Comp_sig_sum = data.table("sig" = unique(Comp$snp))
Comp_sig_sum$N_transeGene = sapply(Comp_sig_sum$sig, function(x) {length(unique(Comp[Comp$snp %in% x, GeneSymbol]))})
Comp_sig_sum$if_transPCO = Comp_sig_sum$sig %in% Us_sig$sig
Comp_sig_sum$trans_eGene = sapply(Comp_sig_sum$sig, function(x) {paste0(unique(Comp[Comp$snp %in% x, GeneSymbol]), collapse = ";")})
Comp_sig_sum$num_included = sapply(strsplit(Comp_sig_sum[, trans_eGene], ";"), function(x) sum(!meta[match(x, meta$gene), Remove], na.rm = TRUE))
Comp_sig_sum = Comp_sig_sum[order(Comp_sig_sum$N_transeGene, decreasing = T), ]

Comp_sig_sum[!Comp_sig_sum$if_transPCO, -'trans_eGene'][1:50, ]
Comp_sig_sum[Comp_sig_sum$if_transPCO, -'trans_eGene']

coexp.module = readRDS("/project2/xuanyao/llw/eQTLGen_DGN_PCO.lambda.01/result/coexp.module.rds")
Us_res = fread("/project2/xuanyao/llw/eQTLGen_DGN_PCO.lambda.01/FDR/signals.chr.module.perm10.txt", header = FALSE, col.names = c("module_snp", "p", "q"))
Us_res$snp = sapply(Us_res$module_snp, function(x) paste0(strsplit(x, ":")[[1]][-1], collapse = ":"))
Us_res$module = sapply(Us_res$module_snp, function(x) strsplit(x, ":")[[1]][1])

snp_of_int = '4:103434253'
trans_target_gene = strsplit(Comp_sig_sum[Comp_sig_sum$sig %in% snp_of_int, trans_eGene], ";")[[1]]

sort(table(coexp.module$moduleLabels[trans_target_gene]), decreasing = TRUE)

Us_res[Us_res$snp %in% snp_of_int, ]
