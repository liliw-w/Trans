rm(list = ls())
require(data.table)

### READ DATA ###
snp_used_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/eQTLGen.used_snp.meta.txt"
gene_used_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/result/DGN.txt"
eqtlGen_res_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz"
eqtlGen_transPCO_res_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/FDR/signals.qvalue.txt"
coexp_module_file = "/scratch/midway2/liliw1/eQTLGen_PCO.lambda.01/result/coexp.module.rds"

snp_used = fread(snp_used_file, header = TRUE)
gene_used = fread(gene_used_file, header = TRUE); gene_used = gene_used[!gene_used$Remove, ]
eqtlGen_res = fread(eqtlGen_res_file, header = TRUE)
eqtlGen_transPCO_res = fread(eqtlGen_transPCO_res_file, header = TRUE)
module = readRDS(coexp_module_file)$moduleLabels

Nmodule = max(module)

### filter eqtlGen's SNPs and genes
ind_overlap = eqtlGen_res$SNP %in% snp_used$SNP & eqtlGen_res$GeneSymbol %in% gene_used$gene
eqtlGen_res = eqtlGen_res[ind_overlap, ]
eqtlGen_res$module = module[match(eqtlGen_res$GeneSymbol, names(module))]

### eqtlGen_transPCO signals under Bonferroni correction
eqtlGen_transPCO_res$module = sapply(eqtlGen_transPCO_res$snp, function(x) paste0("M", strsplit(strsplit(x, ":", fixed = TRUE)[[1]][1], "module")[[1]][2]) )
eqtlGen_transPCO_res$SNP = sapply(eqtlGen_transPCO_res$snp, function(x) strsplit(x, ":", fixed = TRUE)[[1]][2])
eqtlGen_transPCO_res = eqtlGen_transPCO_res[eqtlGen_transPCO_res$p < 0.05/nrow(snp_used)/Nmodule, ]
eqtlGen_transPCO_res$if_module_signal = TRUE

### Organize result for plot
res = do.call(rbind, lapply(unique(eqtlGen_res$SNP), function(x){
  x_table = table(eqtlGen_res[eqtlGen_res$SNP %in% x, module])
  x_num = x_table[as.character(0:Nmodule)]
  x_num[is.na(x_num)] = 0
  x_num = setNames(x_num, paste0("M", 0:Nmodule))
}))
rownames(res) = unique(eqtlGen_res$SNP)

res = as.data.table(res, keep.rownames=TRUE); colnames(res)[1] = "SNP"
res = cbind("SNP_meta" = snp_used[match(res$SNP, snp_used$SNP), meta], res)
res = melt(res, id.vars = c("SNP_meta", "SNP"), variable.name = "module", variable.factor = FALSE, value.name = "module_gene_num")
res = merge(res, eqtlGen_transPCO_res[, 4:6], by = c("SNP", "module"), all.x = TRUE)
res$if_module_signal[is.na(res$if_module_signal)] = FALSE
res = res[order(res$module_gene_num, decreasing = TRUE), ]

res = readRDS('~/Downloads/tmp_res.rds')

Data = res[res$module != "M1"][1:200, ]
Data = res[res$module_gene_num > 5, ]
sum(!Data$if_module_signal)
str(unique(Data[!Data$if_module_signal, SNP]))

Data = res[res$module_gene_num == 0, ]
sum(Data$if_module_signal)
str(unique(Data[Data$if_module_signal, SNP]))
table(Data[Data$if_module_signal, module])

type = "TIE"
if(type == "power"){
  Data = res[res$module_gene_num > 2, ]
  dim(Data)
  cat(sum(!Data$if_module_signal), "(module, SNP) pairs are missed!", '\n')
  cat(length(unique(Data[!Data$if_module_signal, SNP])), "unique SNPs are missed!", '\n')
  ind_point = !Data$if_module_signal; color_point = "red"
  fig_file = "power.png"
}else if(type == "TIE"){
  Data = res[res$module_gene_num <= 2, ]
  dim(Data)
  cat(sum(Data$if_module_signal), "(module, SNP) pairs are falsely selected!", '\n')
  cat(length(unique(Data[Data$if_module_signal, SNP])), "unique SNPs are falsely selected!", '\n')
  ind_point = Data$if_module_signal; color_point = "green"
  fig_file = "TIE.png"
}
p <- ggplot(Data, aes(SNP, module, fill=module_gene_num, label= module_gene_num)) + geom_tile() +
  geom_text(data = Data[Data$module_gene_num!=0, ], size = 1) +
  geom_point(data = Data[ind_point, ], shape = 1, size = 1.5, color = color_point) +
  scale_fill_gradient(low="white", high="grey3") +
  labs(fill = "#") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(fig_file, p, width = 30, height = 10)
