#######
# Check what zscores of snps, whose pvalues are 0, are like. For eQTLGen.
#######

require(data.table)
file.ex.var.regressed = "result/ex.var.regressed.rds"
file.gene.meta = "result/DGN.txt"
file.coexp.module = "result/coexp.module.rds"
chr = 1
module = 1

meta_used_snp = fread('eQTLGen.used_snp.meta.txt', header = TRUE)
res = fread("FDR/signals.qvalue.txt", header = TRUE, col.names = c('module_snp', 'p', 'q'))
res$module = sapply(res$module_snp, function(x) strsplit(x, ":")[[1]][1])
res$snp_ID = sapply(res$module_snp, function(x) strsplit(x, ":")[[1]][2])
res$snp = meta_used_snp[match(res$snp_ID, meta_used_snp$SNP), meta]
res$snp_chr = sapply(res$snp, function(x) strsplit(x, ":")[[1]][1])
#res$snp = sapply(res$module_snp, function(x) paste0(strsplit(x, ":")[[1]][-1], collapse = ":"))
res = as.data.frame(res[order(res$p), ])

datExpr = readRDS(file.ex.var.regressed)
gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module = readRDS(file.coexp.module)
z.mat = fread(paste0('z/z.module', module, '.chr', chr, '.txt.gz'))
z.mat = as.matrix(z.mat, rownames = TRUE)
colnames(z.mat) = gene.meta[match(colnames(z.mat), gene.meta$GeneNameConv), "gene"]
gene_in_cluster = data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene.meta)
gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]
z.mat_trans = z.mat[, gene_trans]; rm(z.mat)

res0 = res[res$p == 0 & res$module == paste0("module", module) & res$snp_chr == chr, ]
res1 = res[res$p != 0 & res$module == paste0("module", module) & res$snp_chr == chr, ]; res1 = res1[order(res1$p, decreasing = TRUE), ]
z0 = z.mat_trans[match(res0$snp_ID, rownames(z.mat_trans)), ]
z1 = z.mat_trans[match(res1$snp_ID, rownames(z.mat_trans)), ]

## Check summary
sort(table(res[res$p == 0, 'module']), decreasing = TRUE)
sort(table(res[res$p == 0 & res$module == "module1", "snp_chr"]), decreasing = TRUE)
sort(table(meta_used_snp$SNPChr), decreasing = TRUE)

## Plot zscores for each gene
pdf('tmp.pdf', width = 5*3, height = 5)
ord = order(abs(z0[1, ]), decreasing = T)
par(mfrow = c(1, 3))
matplot(t(z0[, ord]), type = "p", pch = 1, col = 2, xlab = 'Gene', ylab = "z-score", ylim = c(-20, 20))
matpoints(t(z1[sample(nrow(z1), 10), ord]), type = "p", pch = 2, col = 1)
hist(z0[sample(nrow(z0), 1), ], xlim = c(-10, 10), ylim = c(0, 80), breaks = seq(-100, 100, 0.5), xlab = "z of (SNP, gene) with 0 p", ylab = "#Genes in module", main = paste0("module", module))
hist(z1[sample(nrow(z1), 1), ], xlim = c(-10, 10), ylim = c(0, 80), breaks = seq(-100, 100, 0.5), xlab = "z of (SNP, gene) with non-0 p", ylab = "#Genes in module", main = paste0("module", module))
dev.off()

## Plot zscores for each SNP
pdf('tmp.pdf', width = 5, height = 5)
matplot(z0[, ], type = "p", pch = 1, col = 2, xlab = 'SNP', ylab = "z-score", ylim = c(-20, 20))
dev.off()

## plot #genes with large z-scores for each SNP
res_SNP = NULL
for(thre_p_z in c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)){
  thre_z2 = qchisq(1-thre_p_z, 1)
  num_large_z = rowSums(z0^2 > thre_z2)
  res_SNP = rbind(res_SNP, data.table("num_large_z" = sort(num_large_z, decreasing = TRUE),
                                      "thre_z" = rep(thre_p_z, nrow(z0)),
                                      "SNP" = 1:nrow(z0)))
}
fig = ggplot(res_SNP, aes(x = SNP, y = num_large_z, color = factor(thre_z))) +
  geom_point(size = 0.9) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(vars(thre_z))
ggsave("large.z.png", fig)


## plot #SNPs on chr's included v.s. #SNPs on chr's with 0 pvalues corresponding to the module
png('tmp.png')
plot(as.numeric(table(meta_used_snp$SNPChr)[as.character(1:22)]),
     as.numeric(table(res[res$p == 0 & res$module == paste0("module", module), "snp_chr"])[as.character(1:22)]),
     xlab = "#snps on Chromosome in eQTLGen",
     ylab = "#snps corresponding to module with 0 pvalues")
dev.off()

