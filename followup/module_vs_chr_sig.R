rm(list = ls())
require(data.table)
require(ggplot2)
require(tidyverse)

## Input files
file_signals = "/scratch/midway2/liliw1/MODULES/MSigDB/FDR/signals.chr.module.perm10.fdr10.txt"
file_indep_signals = "/scratch/midway2/liliw1/MODULES/MSigDB/postanalysis/indep.signals.chr.module.perm10.fdr10.txt"
file_coexp_module = "/scratch/midway2/liliw1/MODULES/MSigDB/result/coexp.module.rds"
if_module_anno = TRUE
plot_path = "./plot/"
Nchr = 22

## read files
signals = fread(file_signals, header = FALSE)
indep_signals = fread(file_indep_signals, header = FALSE)
coexp_module = readRDS(file_coexp_module)

## summarize data
Nmodule = max(coexp_module$moduleLabels)
module_size = as.numeric(table(coexp_module$moduleLabels))
signals$module = sapply(signals$V1, function(x) strsplit(x, ":", fixed = TRUE)[[1]][1])
signals$module_id = as.numeric(sapply(signals$module, function(x) strsplit(x, "module", fixed = TRUE)[[1]][2]))
signals$chr = as.numeric(sapply(signals$V1, function(x) strsplit(x, ":", fixed = TRUE)[[1]][2]))

### module annotaions
if(if_module_anno){
  ## data files
  file_gene_set = "/scratch/midway2/liliw1/MODULES/MSigDB/h.all.v7.4.symbols.gmt"
  file.ex.var.regressed = "~/xuanyao_llw/DGN_no_filter_on_mappability/result/ex.var.regressed.rds"

  ## read files
  gene_set = qusage::read.gmt(file_gene_set)
  datExpr = readRDS(file.ex.var.regressed)

  ## gene sets with genes included in data
  genes_data = colnames(datExpr)
  gene_set_in_data = lapply(gene_set, function(x) x[x %in% genes_data])
  module_size = sort(sapply(gene_set_in_data, function(x) length(x)), decreasing = TRUE)
  module_anno = names(module_size)
}else{module_anno = NULL}


## data frame for plotting
df_signals = signals %>% group_by(module_id, chr) %>% summarise(num_sig = n())
df_signals$module_id = factor(df_signals$module_id, levels = 1:Nmodule,
                              labels = paste0("M", 1:Nmodule, "(#", module_size, ")", "-", module_anno) )
df_signals$chr = factor(df_signals$chr, levels = 1:Nchr, labels = paste0("Chr", 1:Nchr))


### tile plot
fig_tile <- ggplot(df_signals, aes(chr, module_id)) +
  geom_tile(aes(fill = num_sig)) +
  scale_fill_gradient(low="white", high="grey2") +
  geom_text(data = df_signals[df_signals$num_sig>0, ], aes(label = num_sig), color = "tomato3", size = 2) +
  labs(x = "Chromosome", y = "Module", fill = "#signals") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title=element_text(size=5))
ggsave(paste0("num_signals_module_vs_chr_", ".png"), fig_tile, path = plot_path)
system(paste0("bash ~/imgcat ", plot_path, "num_signals_module_vs_chr_", ".png"))


