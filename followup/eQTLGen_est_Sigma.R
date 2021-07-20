rm(list = ls())
require(data.table)

meta=fread('~/scratch/eQTLGen_PCO.lambda.01/result/DGN.txt')
module_info=readRDS('~/scratch/eQTLGen_PCO.lambda.01/result/coexp.module.rds')$moduleLabels
names(module_info) = meta[match(names(module_info), meta$gene), GeneNameConv]
Nmodule = max(module_info)

## Extract (trans) summary stat for each modules
a = fread('/project2/xuanyao/data/eQTLGen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz')
for(module in 1:Nmodule){
  module1 = names(module_info)[module_info==module]
  a_module1 = a[a$Gene %in% module1, ]

  saveRDS(a_module1, paste0("sumStat.trans.module", module, ".rds"))
  cat("module:", module, "\n")
}


## Extract (cis) summary stat for each modules
a = fread('/project2/xuanyao/data/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')

for(module in 10:Nmodule){
  module1 = names(module_info)[module_info==module]
  a_module1 = a[a$Gene %in% module1, ]

  saveRDS(a_module1, paste0("sumStat.cis.module", module, ".rds"))
  cat("module:", module, "\n")
}


## Combine trans and cis summary stat for each module
## And extract the common SNPs, i.e. SNPs with summary stat for all genes within each module.
for(module in 1:Nmodule){
  a_module1 = readRDS(paste0("sumStat.cis.module", module, ".rds"))
  b_module1 = readRDS(paste0("sumStat.trans.module", module, ".rds"))
  ab_module1 = rbind(a_module1, b_module1)
  ab_module1$SNPmeta = paste(ab_module1$SNPChr, ab_module1$SNPPos, sep = ":")

  ab_z = dcast(ab_module1, SNP+SNPmeta~Gene, fun.aggregate = max,
               value.var = 'Zscore',
               drop = TRUE, fill = NA)
  ab_z_na = is.na(ab_z)
  rem_ind = apply(ab_z_na[, -c(1:2)], 1, sum )
  ab_z_remain = ab_z[rem_ind == 0, ]
  saveRDS(ab_z_remain, paste0("sumStat.all.module", module, ".rds"))

  cat("module:", module, "\n")
}


## #Indep null SNPs (and #Uniq null SNPs) for each module under various z-score threshold for being null
res_nullSNP = NULL
for(module in 1:Nmodule){
  ab_z_remain = readRDS(paste0("sumStat.all.module", module, ".rds"))
  NSNP = nrow(ab_z_remain); module_size = ncol(ab_z_remain)-2

  for(thre_p_z in c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)){
    thre_z2 = qchisq(1-thre_p_z, 1)
    altSNP_ind = apply((ab_z_remain[, -c(1:2)])^2 > thre_z2, 1, sum) > 0

    file_uniq = paste0("eQTLGen_est_Sigma/unique.SNP.module", module, ".", thre_p_z, ".txt")
    file_indep = paste0("eQTLGen_est_Sigma/indep.SNP.module", module, ".", thre_p_z, ".txt")
    fwrite(ab_z_remain[!altSNP_ind, "SNPmeta"], file_uniq,
           quote = FALSE, sep = "\t", col.names = FALSE)
    command = paste0("bash eQTLGen_est_Sigma/indep.SNP.sh ", file_uniq, " ", file_indep)
    system(command)
    num_nullSNP_indep = as.numeric(system(paste0("cat ", file_indep, " | wc -l"), intern = TRUE))

    res_nullSNP = rbind(res_nullSNP, data.table("thre_z" = thre_p_z,
                                                "num_nullSNP_uniq" = NSNP - sum(altSNP_ind),
                                                "num_nullSNP_indep" = num_nullSNP_indep,
                                                "module" = module,
                                                "module_size" = module_size))
  }
}
saveRDS(res_nullSNP, "num_nullSNP.rds")


## Plot #Indep null SNPs (and #Uniq null SNPs) for each module under various z-score threshold for being null
res_nullSNP = readRDS('num_nullSNP.rds')
require(ggplot2)
fig = ggplot(data = res_nullSNP, aes(x = factor(module)) ) +
  geom_bar(aes(y = -module_size), stat="identity") +
  geom_point(aes(y = num_nullSNP_uniq, color = factor(thre_z), shape = factor("Uniq null SNP")), alpha = 0.9) +
  geom_point(aes(y = num_nullSNP_indep, color = factor(thre_z), shape = factor("Indep null SNP")) ) +
  geom_text(data = res_nullSNP[!duplicated(res_nullSNP$module), ], aes(y = -module_size*10-100, label = module_size)) +
  scale_shape_manual(values = c(19, 3)) +
  #ylim(0, 10000) +
  labs(title = "#null SNP", x = "module", y = "num_nullSNP", color="thre_z", shape = "SNP") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave('num_nullSNP.png', fig, height = 10, width = 15)



## Check how Sigma from (1) independent null eQTLGen z-scores , (2) DGN expression data, affect PCO pvalues
require(data.table)
### READ DATA ###
eqtlGen_transPCO_res_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/FDR/signals.qvalue.txt"
coexp_module_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/result/coexp.module.rds"
file_eQTLGen_z = "/scratch/midway2/liliw1/eQTLGen_est_Sigma/sumStat.all.module1.rds"
file_indep_null_snp = "/scratch/midway2/liliw1/eQTLGen_est_Sigma/indep.SNP.module1.0.001.txt"
gene_meta_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/result/DGN.txt"

eqtlGen_transPCO_res = fread(eqtlGen_transPCO_res_file, header = TRUE)
module_info = readRDS(coexp_module_file)$moduleLabels; Nmodule = max(module_info)
eQTLGen_z = readRDS(file_eQTLGen_z)
indep_null_snp = fread(file_indep_null_snp, header = FALSE)$V1
gene.meta = fread(gene_meta_file, header = TRUE)

### Estimate the correlation matrix of zscore using extracted independent null z's, for the whole module ###
indep_null_z = eQTLGen_z[eQTLGen_z$SNPmeta %in% indep_null_snp, -c(1:2)]
Sigma_null_z = cor(indep_null_z)

### Apply PCO using this estimate Sigma ###
params1 = "/project2/xuanyao/llw/TCGA_PCO.lambda.01/script/"
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))

file.ex.var.regressed = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/result/ex.var.regressed.rds"
datExpr = readRDS(file.ex.var.regressed)

module = 1
gene_in_cluster = data.frame("gene" = names(module_info[module_info == module]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene.meta)
for (chr in 1:22) {
  file.z = paste0("/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/z/z.module", module, ".chr", chr, ".txt.gz")
  z.mat = fread(file.z)
  z.mat = as.matrix(z.mat, rownames = TRUE)
  gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "GeneNameConv"]
  z.mat_trans = z.mat[, gene_trans]

  ### Use Sigma from eQTLGen independent null zscores ###
  Sigma = Sigma_null_z[match(gene_trans, rownames(Sigma_null_z)), match(gene_trans, colnames(Sigma_null_z))]
  SigmaO = ModifiedSigmaOEstimate(Sigma)
  p.all = ModifiedPCOMerged(Z.mat=z.mat_trans,
                            Sigma=Sigma, SigmaO=SigmaO)

  ### Use Sigma from DGN expression data ###
  Sigma_DGN = cor(datExpr[, gene.meta[match(gene_trans, gene.meta$GeneNameConv), gene]])
  SigmaO_DGN = ModifiedSigmaOEstimate(Sigma_DGN)
  p.all_DGN = ModifiedPCOMerged(Z.mat=z.mat_trans,
                                Sigma=Sigma_DGN, SigmaO=SigmaO_DGN)

  ### Organize results ###
  res = data.table(data.frame("p_nullz" = p.all,
                              "p_DGN" = p.all_DGN, check.rows = TRUE),
                   keep.rownames = TRUE)
  res = res[order(res$p_DGN), ]
  saveRDS(res, paste0('p.module', module, ".chr", chr, '.Sigma.rds'))

  cat('module', module, "chr", chr, "has done!", "\n")
}


## Plot PCO pvalues using the two Sigma's
require(data.table)
require(ggplot2)
library(gridExtra)
snp_used_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/eQTLGen.used_snp.meta.txt"
coexp_module_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/result/coexp.module.rds"
snp_used = fread(snp_used_file, header = TRUE)
module_info = readRDS(coexp_module_file)$moduleLabels; Nmodule = max(module_info)

module = 1
res_all = NULL
for (chr in 1:22){
  res = readRDS(paste0('p.module', module, ".chr", chr, '.Sigma.rds'))
  res$x = 1:nrow(res)
  res$chr = chr
  res_all = rbind(res_all, res)
}

### Plot and compare the PCO pvalues for this module using the two Sigma's ###
fig_p = ggplot(data = res_all, aes(x = x)) +
  geom_point(aes(y = log10(p_nullz), color = factor("nullz")), size = 0.5 ) +
  geom_point(aes(y = log10(p_DGN), color = factor("DGNz")), size = 0.5 ) +
  geom_hline(yintercept = log10(0.05/nrow(snp_used)/Nmodule), color = "grey", linetype = "dashed") +
  facet_wrap(vars(chr), scales = "free_x", labeller = "label_both") +
  labs(title = paste0('Effect of Sigma on PCO pvalue, module', module), x = "SNP on a chr", y = "Log10(p) of a module", color="z") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste0('effect.of.Sigma.on.PCO.pvalue.module', module, '.png'), fig_p,
       height = 10, width = 10)
# system(paste0('bash ~/imgcat effect.of.Sigma.on.PCO.pvalue.module', module, '.png'))

cat("Under Bonferroni correction, pvalue threshold is ", 0.05/nrow(snp_used)/Nmodule, "\n", "#Signals for this module using two Sigma's are: ", "\n")
print(apply(res_all[, c("p_nullz", "p_DGN")], 2, function(x) sum(x < 0.05/nrow(snp_used)/Nmodule)))


### QQ-plot of the PCO pvalues for this module using the two Sigma's ###
qqplot.hist <- function(input, title){
  observed <- sort(input)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected+1))))
  df = data.frame(x = lexp, y = lobs, yy = observed)

  res = list()
  res[[1]] = ggplot(df, aes(x=x, y=y)) + geom_point(size = 0.2) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    labs(title = title, x = "Expected (-logP)", y = "Observed (-logP)") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  res[[2]] = ggplot(df, aes(yy)) + geom_histogram(breaks = seq(0, 1, 0.01)) +
    labs(title = title, x = "Observed (P)") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(res)
}

old_p = readRDS('../eQTLGen_PCO.lambda.01/p.obs.all.rds')
# qqplot.hist(old_p[old_p$module == paste0("module", module), "p"], "old_DGNz")

fig_qq = rbind(qqplot.hist(res_all$p_nullz, "nullz"), qqplot.hist(res_all$p_DGN, "DGNz") )
ggsave(paste0('qqplot.effect.of.Sigma.on.PCO.pvalue.module', module, '.png'),
       marrangeGrob(fig_qq, nrow=length(fig_qq)/2, ncol=2, top = NULL),
       height = length(fig_qq)/2*(6/2), width = 6)
# system(paste0('bash ~/imgcat qqplot.effect.of.Sigma.on.PCO.pvalue.module', module, '.png'))

#### eqtlGen_transPCO signals under Bonferroni correction
eqtlGen_transPCO_res$module = sapply(eqtlGen_transPCO_res$snp, function(x) strsplit(x, ":", fixed = TRUE)[[1]][1] )
eqtlGen_transPCO_res$SNP = sapply(eqtlGen_transPCO_res$snp, function(x) strsplit(x, ":", fixed = TRUE)[[1]][2])
eqtlGen_transPCO_res$q = eqtlGen_transPCO_res$p*nrow(snp_used)*Nmodule
eqtlGen_transPCO_res = eqtlGen_transPCO_res[eqtlGen_transPCO_res$q < 0.05, ]

