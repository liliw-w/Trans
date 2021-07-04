#######
# Check what zscores of snps, whose pvalues are 0, are like. For eQTLGen.
#######

rm(list = ls())
require(data.table)
require(ggplot2)
library(gridExtra)
require(tidyverse)
file.gene.meta = "result/gene.meta.txt" # "result/DGN.txt" # "result/gene.meta.txt"
file.coexp.module = "result/coexp.module.rds"
file.signal = "FDR/signals.chr.module.perm10.txt" # "FDR/signals.qvalue.txt" # "FDR/signals.chr.module.perm10.txt"
is_eQTGen = FALSE
module = 1
chr = 1

## Data preparation ##
p_all = data.table(data.frame("p" = readRDS(paste0("p/p.module", module, ".chr", chr, ".rds"))), keep.rownames = TRUE)
signal = as.data.frame(fread(file.signal, header = TRUE, col.names = c('module_snp', 'p', 'q')))

res = as.data.frame(data.table("module_snp" = paste(paste0("module", module), p_all$rn, sep = ":"),
                               "module" = rep(paste0("module", module), nrow(p_all)),
                               "snp_ID" = p_all$rn,
                               "p" = p_all$p,
                               #"snp" = meta_used_snp[match(p_all$rn, meta_used_snp$SNP), meta],
                               "snp_chr" = rep(chr, nrow(p_all))))
res$if_sig = res$module_snp %in% signal$module_snp
if(!is_eQTGen){
  ind_large_p = which(res$if_sig != TRUE & res$p < 1e-4)
  ind_large_p = c(ind_large_p, sample(which(res$if_sig != TRUE & res$p > 1e-4 & res$p < 1e-3), 200))
  res = rbind(res[res$if_sig, ], res[ind_large_p, ])
}
res = res[order(res$p), ]


gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module = readRDS(file.coexp.module)
z.mat = fread(paste0('z/z.module', module, '.chr', chr, '.txt.gz'))
z.mat = as.matrix(z.mat, rownames = TRUE)
if(is_eQTGen){colnames(z.mat) = gene.meta[match(colnames(z.mat), gene.meta$GeneNameConv), "gene"]}
gene_in_cluster = data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene.meta)
gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]
z.mat_trans = z.mat[rownames(z.mat) %in% res$snp_ID, gene_trans]; rm(z.mat)


## Plot zscores for each SNP ##
DT_z.mat_trans = data.table(z.mat_trans, keep.rownames = TRUE)
DT_z.mat_trans$p = res[match(DT_z.mat_trans$rn, res$snp_ID), 'p']
DT_z.mat_trans$if_sig = res[match(DT_z.mat_trans$rn, res$snp_ID), 'if_sig']
DT_z.mat_trans$snp_x = nrow(DT_z.mat_trans)+1-rank(DT_z.mat_trans$p, ties.method = "first")
df_zmat = melt(DT_z.mat_trans, id.vars = c("rn", "p", 'if_sig', "snp_x"), variable.name = "gene", value.name = "z", variable.factor = FALSE)

fig_snp = ggplot(df_zmat[, ], aes(x = snp_x, y = z, color = factor(round(-log10(p))), shape = factor(ifelse(if_sig, "sig", "non-sig")))) +
  geom_point(size = 0.7) +
  ylim(-20, 20) +
  labs(title = "zscores", x = "SNP", y = "z-score", color="-log(p)", shape = "if_sig") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

fig_abs_z = ggplot(df_zmat[, ], aes(x = -log10(p), y = abs(z), color = factor(round(-log10(p))), shape = factor(ifelse(if_sig, "sig", "non-sig")))) +
  geom_point(size = 0.7) +
  ylim(-20, 20) +
  labs(title = "abs zscores", x = "-log10(p)", y = "abs(z-score)", color="-log(p)", shape = "if_sig") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Plot genes with large zscores for each SNP under various cutoff ##
res_SNP_numGene = NULL
for(thre_p_z in c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)){
  thre_z2 = qchisq(1-thre_p_z, 1)
  res_SNP_numGene = rbind(res_SNP_numGene, df_zmat %>% group_by(rn) %>% summarise("p" = unique(p),
                                                                                  "if_sig" = unique(if_sig),
                                                                                  "snp_x" = unique(snp_x),
                                                                                  "num_large_z" = sum(z^2 > thre_z2),
                                                                                  "thre_z" = thre_p_z) %>%
                            group_by("p_interval" = round(-log10(p))) %>%
                            mutate("group_order" = (min(snp_x):max(snp_x))[rank(num_large_z, ties.method = "first")]))
}

fig_numGene = ggplot(res_SNP_numGene, aes(x = group_order, y = num_large_z, color = factor(p_interval), shape = factor(ifelse(if_sig, "sig", "non-sig")))) +
  geom_point(size = 0.7) +
  facet_wrap(vars(thre_z), labeller = "label_both") +
  labs(title = "large zscores", x = "SNP", y = "num_large_z", color="-log(p)", shape = "if_sig") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(paste0("check_z_of_0p.module", module, ".chr", chr, ".png"),
       marrangeGrob(rbind(list(fig_snp), list(fig_abs_z), list(fig_numGene)), nrow=3, ncol=1, top = NULL),
       height = 18, width = 6)


run_more = FALSE
if(run_more & is_eQTGen){
  meta_used_snp = fread('eQTLGen.used_snp.meta.txt', header = TRUE)

  signal = as.data.frame(fread(file.signal, header = TRUE, col.names = c('module_snp', 'p', 'q')))
  signal$module = sapply(signal$module_snp, function(x) strsplit(x, ":")[[1]][1])
  signal$snp_ID = sapply(signal$module_snp, function(x) strsplit(x, ":")[[1]][2])
  signal$snp = meta_used_snp[match(signal$snp_ID, meta_used_snp$SNP), meta]
  signal$snp_chr = sapply(signal$snp, function(x) strsplit(x, ":")[[1]][1])
  #signal$snp = sapply(signal$module_snp, function(x) paste0(strsplit(x, ":")[[1]][-1], collapse = ":"))


  ## Plot zscores for each gene
  res0 = res[res$p == 0, ]
  res1 = res[res$p > 1e-7, ]; res1 = res1[order(res1$p, decreasing = TRUE), ]
  z0 = z.mat_trans[match(res0$snp_ID, rownames(z.mat_trans)), ]
  z1 = z.mat_trans[match(res1$snp_ID, rownames(z.mat_trans)), ]

  pdf('z_xgene.pdf', width = 5*3, height = 5)
  ord = order(abs(z0[1, ]), decreasing = T)
  par(mfrow = c(1, 3))
  matplot(t(z0[, ord]), type = "p", pch = 1, col = 2, xlab = 'Gene', ylab = "z-score", ylim = c(-20, 20))
  matpoints(t(z1[, ord]), type = "p", pch = 2, col = 1)
  hist(z0[sample(nrow(z0), 1), ], xlim = c(-10, 10), ylim = c(0, 80), breaks = seq(-100, 100, 0.5), xlab = "z of (SNP, gene) with 0 p", ylab = "#Genes in module", main = paste0("module", module))
  hist(z1[sample(nrow(z1), 1), ], xlim = c(-10, 10), ylim = c(0, 80), breaks = seq(-100, 100, 0.5), xlab = "z of (SNP, gene) with non-0 p", ylab = "#Genes in module", main = paste0("module", module))
  dev.off()


  ## plot #SNPs on chr's included v.s. #SNPs on chr's with 0 pvalues corresponding to the module
  png('tmp.png')
  plot(as.numeric(table(meta_used_snp$SNPChr)[as.character(1:22)]),
       as.numeric(table(res[res$p == 0 & res$module == paste0("module", module), "snp_chr"])[as.character(1:22)]),
       xlab = "#snps on Chromosome in eQTLGen",
       ylab = "#snps corresponding to module with 0 pvalues")
  dev.off()


  ## Check summary
  module = 1
  chr = 1

  sort(table(signal$module), decreasing = TRUE)
  module = 11
  sort(table(signal[signal$module == paste0("module", module), "snp_chr"]), decreasing = TRUE)
  sort(table(signal[signal$module == paste0("module", module) & signal$p < 1e-7, "snp_chr"]), decreasing = TRUE)
  sort(table(meta_used_snp$SNPChr), decreasing = TRUE)

  sort(table(signal[signal$p == 0, 'module']), decreasing = TRUE)
  sort(table(signal[(signal$p == 0 & signal$module == paste0("module", module)), "snp_chr"]), decreasing = TRUE)

}

