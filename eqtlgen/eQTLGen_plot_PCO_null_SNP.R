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
  geom_point(aes(y = -log10(p_nullz), color = factor("nullz")), size = 0.5 ) +
  geom_point(aes(y = -log10(p_DGN), color = factor("DGNz")), size = 0.5 ) +
  #geom_hline(yintercept = log10(0.05/nrow(snp_used)/Nmodule), color = "grey", linetype = "dashed") +
  #facet_wrap(vars(chr), scales = "free_x", labeller = "label_both") +
  labs(title = paste0('Effect of Sigma on PCO pvalue, module', module), x = "SNP on a chr", y = "-Log10(p) of a module", color="z") +
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
