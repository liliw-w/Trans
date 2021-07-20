rm(list = ls())
require(data.table)
file_module = '/project2/xuanyao/llw/DGN_PCO.lambda.01/result/coexp.module.rds'

module_info = readRDS(file_module)$moduleLabels; Nmodule = max(module_info)
thre_p_z = 1e-4
thre_z2 = qchisq(1-thre_p_z, 1)

## Indep null SNPs (and #Uniq null SNPs) for each module under various z-score threshold for being null
for(module in 1:Nmodule){
  for(chr in 1:22){
    #### read z-scores ####
    file_zmat = paste0('/project2/xuanyao/llw/DGN_PCO.lambda.01/z/z.module', module, '.chr', chr, '.txt.gz')
    zmat = fread(file_zmat, header = TRUE)

    #### extract SNPs with small z-scores (null SNPs) ####
    file_uniq = paste0("unique.SNP/unique.SNP.module", module, '.chr', chr, ".", thre_p_z, ".txt")
    altSNP_ind = apply((zmat[, -c(1)])^2 > thre_z2, 1, sum) > 0
    fwrite(zmat[!altSNP_ind, "snp"], file_uniq,
           quote = FALSE, sep = "\t", col.names = FALSE)

    #### extract independent null SNPs ####
    file_indep = paste0("indep.SNP/indep.SNP.module", module, '.chr', chr, ".", thre_p_z, ".txt")
    command = paste0("bash indep.SNP.sh ", file_uniq, " ", file_indep, " ", chr)
    system(command)

    #### write out zscores of these independent null SNPs ####
    file_zmat_indep_null = paste0("z.indep.null/z.indep.null.module", module, '.chr', chr, ".", thre_p_z, ".txt.gz")
    snp_indep_null = fread(file_indep, header = FALSE)
    zmat_indep_null = zmat[zmat$snp %in% snp_indep_null$V1, ]
    fwrite(zmat_indep_null, file_zmat_indep_null,
           quote = FALSE, sep = "\t", col.names = TRUE)

    cat("Chromosome:", chr, "\n")
  }
  #### bind all indep null SNPs ####
  file_zmat_indep_null = paste0("z.indep.null/z.indep.null.module", module, '.chr', 1:22, ".", thre_p_z, ".txt.gz")
  zmat_indep_null = rbindlist(lapply(file_zmat_indep_null, function(x) fread(x, header = TRUE)),
                              use.names = TRUE)
  saveRDS(zmat_indep_null, paste0("z.indep.null/combined.z.indep.null.module", module, ".", thre_p_z, ".rds"))

  #### calculate sigma using these independent null z-scores ####
  Sigma_null_z = cor(zmat_indep_null[, -'snp'])
  saveRDS(Sigma_null_z, paste0("Sigma.null.z/Sigma.null.z.module", module, ".", thre_p_z, ".rds"))
}


## Check how Sigma from (1) independent null z-scores , (2) DGN expression data, affect PCO pvalues
command = 'bash Snakefile_DGN_est_Sigma.sh'
system(command)


## plot distribution of pvalues using Sigma_nullz and Sigma_exp
require(data.table)
require(ggplot2)
library(gridExtra)
qqplot.hist <- function(input, title){
  observed <- sort(input)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected+1))))
  df = data.frame(x = lexp, y = lobs, yy = observed)

  res = list()
  res[[1]] = ggplot(df, aes(x=x, y=y)) + geom_point(size = 0.1) +
    geom_abline(slope = 1, intercept = 0, color="red") +
    labs(title = title, x = "Expected (-logP)", y = "Observed (-logP)") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  res[[2]] = ggplot(df, aes(yy)) + geom_histogram(breaks = seq(0, 1, 0.01)) +
    labs(title = title, x = "Observed (P)") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(res)
}


module = 1
res_all = NULL
for (chr in 1:22){
  file_p_nullz = paste0('p.Sigma.null.z/p.module', module, ".chr", chr, '.rds')
  file_p_DGN = paste0('/project2/xuanyao/llw/DGN_PCO.lambda.01/p/p.module', module, ".chr", chr, '.rds')

  p.all_nullz = readRDS(file_p_nullz)
  p.all_DGN = readRDS(file_p_DGN)
  res = data.table(data.frame("p_nullz" = p.all_nullz,
                              "p_DGN" = p.all_DGN, check.rows = TRUE),
                   keep.rownames = TRUE)

  res = res[order(res$p_DGN), ]
  res$x = 1:nrow(res)
  res$chr = chr
  res_all = rbind(res_all, res)
}

### plot 1: pvalue_sigma_nullz v.s. pvalue_sigma_exp on log10 scale ###
fig_p = ggplot(data = res_all, aes(x = x)) +
  geom_point(aes(y = log10(p_nullz), color = factor("nullz")), size = 0.1 ) +
  geom_point(aes(y = log10(p_DGN), color = factor("DGNz")), size = 0.1 ) +
  #geom_hline(yintercept = log10(0.05/nrow(snp_used)/Nmodule), color = "grey", linetype = "dashed") +
  facet_wrap(vars(chr), scales = "free_x", labeller = "label_both") +
  labs(title = paste0('Effect of Sigma on PCO pvalue, module', module), x = "SNP on a chr", y = "Log10(p) of a module", color="z") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste0('effect.of.Sigma.on.PCO.pvalue.module', module, '.log10.png'), fig_p,
       height = 10, width = 10)
system(paste0('bash ~/imgcat effect.of.Sigma.on.PCO.pvalue.module', module, '.log10.png'))


### plot 2: pvalue_sigma_nullz v.s. pvalue_sigma_exp on orignal scale ###
fig_p = ggplot(data = res_all, aes(x = x)) +
  geom_point(aes(y = (p_nullz), color = factor("nullz")), size = 0.1 ) +
  geom_point(aes(y = (p_DGN), color = factor("DGNz")), size = 0.1 ) +
  #geom_hline(yintercept = log10(0.05/nrow(snp_used)/Nmodule), color = "grey", linetype = "dashed") +
  facet_wrap(vars(chr), scales = "free_x", labeller = "label_both") +
  labs(title = paste0('Effect of Sigma on PCO pvalue, module', module), x = "SNP on a chr", y = "p of a module", color="z") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste0('effect.of.Sigma.on.PCO.pvalue.module', module, '.origp.png'), fig_p,
       height = 10, width = 10)
system(paste0('bash ~/imgcat effect.of.Sigma.on.PCO.pvalue.module', module, '.origp.png'))


### plot 3: QQplot of pvalue_sigma_nullz and pvalue_sigma_exp ###
fig_qq = rbind(qqplot.hist(res_all$p_nullz, "nullz"), qqplot.hist(res_all$p_DGN, "DGNz") )
ggsave(paste0('qqplot.effect.of.Sigma.on.PCO.pvalue.module', module, '.png'),
       marrangeGrob(fig_qq, nrow=length(fig_qq)/2, ncol=2, top = NULL),
       height = length(fig_qq)/2*(6/2), width = 6)
system(paste0('bash ~/imgcat qqplot.effect.of.Sigma.on.PCO.pvalue.module', module, '.png'))
