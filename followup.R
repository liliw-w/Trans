######## Find all cis-eQTLs for every eGene ###########
require(data.table)

file_input = '~/xuanyao_llw/DGN_PCO.lambda.01/DGN_eQTL.txt.gz'
file_eGene = '~/xuanyao_llw/DGN_PCO.lambda.01/DGN_eQTL.txt.gz'
fdr = 0.05

eGene = fread(file_input, header=TRUE)
eGene = eGene[eGene$qval < fdr, ]

pt = mean(c(max(eGene$bpval), min(eGene$bpval)))
eGene$pval_nominal_threshold = signif(qbeta(pt,
                                            eGene$shape1, eGene$shape2, ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

fwrite(eGene, file_eGene, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


file_eGene = '~/xuanyao_llw/DGN_PCO.lambda.01/DGN_eQTL.txt.gz'
file_all = 'nom_20phenoPC_allChr.txt.gz'
file_eQTL = '~/xuanyao_llw/DGN_PCO.lambda.01/DGN_eQTL_signif_variant_gene_pairs.txt.gz'
fdr = 0.05

eGene = fread(file_eGene, header = TRUE)
assoc_all = fread(file_all, header = TRUE, sep=" ")

eGene = eGene[eGene$qval < fdr, ]
sig_egene = unique(eGene$pid)
assoc_all = assoc_all[assoc_all$pid %in% sig_egene, ]

assoc_all$pval_nominal_threshold = eGene[match(assoc_all$pid, eGene$pid), pval_nominal_threshold]
assoc_all$if_QTL = assoc_all$npval <= assoc_all$pval_nominal_threshold

fwrite(assoc_all[assoc_all$if_QTL, ], file_eQTL, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


a = assoc_all[2000:5000, ]
sum(a$pid %in% sig_egene)
a = a[a$pid %in% sig_egene, ]


library(ggplot2)
library(gridExtra)
library(data.table)


Nmodule = 146; plot.name = "plots/Muscle_Skeleta.p.funcExplorer.png"
file.p = as.character(outer(1:Nmodule, 1:22, FUN = function(x, y) paste0("p/p.module", x, ".chr", y, ".rds")))
p.obs = rbindlist(lapply(file.p, function(x)
{tmp_y=readRDS(x);
print(x);
tmp_y = tmp_y[tmp_y<1e-6];
if(length(tmp_y > 0)){as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)}
}))
p.obs$V2[p.obs$V2==0] = 1e-20


Nmodule = 32; plot.name = "plots/Muscle_Skeleta.p.min20_deep4.png"
file.p = as.character(outer(1:Nmodule, 1:22, FUN = function(x, y) paste0("p_min20_deep4/p.module", x, ".chr", y, ".rds")))
p.obs2 = rbindlist(lapply(file.p, function(x)
{tmp_y=readRDS(x);
tmp_y = tmp_y[tmp_y<1e-6];
if(length(tmp_y > 0)){as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)}
}))
p.obs2$V2[p.obs2$V2==0] = 1e-20


Nmodule = 18; plot.name = "plots/Muscle_Skeleta.p.min20_deep2.png"
file.p = as.character(outer(1:Nmodule, 1:22, FUN = function(x, y) paste0("~/xuanyao_llw/GTEx_v8/Muscle_Skeletal/p/p.module", x, ".chr", y, ".rds")))
p.obs3 = rbindlist(lapply(file.p, function(x)
{tmp_y=readRDS(x);
tmp_y = tmp_y[tmp_y<1e-6];
if(length(tmp_y > 0)){as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)}
}))
p.obs3$V2[p.obs3$V2==0] = 1e-20


fig1 = ggplot(old$'p.min10_deep4', aes(x=-log10(V2))) + geom_histogram(binwidth=0.5) + coord_cartesian(xlim=c(0,20), ylim=c(0, 300)) + ggtitle("Muscle_Skeleta.p.min10_deep4") + xlab("-log10(p)") + theme(text=element_text(size=5))
fig2 = ggplot(old$'p.min20_deep4', aes(x=-log10(V2))) + geom_histogram(binwidth=0.5) + coord_cartesian(xlim=c(0,20), ylim=c(0, 300)) + ggtitle("Muscle_Skeleta.p.min20_deep4") + xlab("-log10(p)") + theme(text=element_text(size=5))
fig3 = ggplot(old$'p.min20_deep2', aes(x=-log10(V2))) + geom_histogram(binwidth=0.5) + coord_cartesian(xlim=c(0,20), ylim=c(0, 300)) + ggtitle("Muscle_Skeleta.p.min20_deep2") + xlab("-log10(p)") + theme(text=element_text(size=5))
fig4 = ggplot(old$'p.funcExplorer', aes(x=-log10(V2))) + geom_histogram(binwidth=0.5) + coord_cartesian(xlim=c(0,20), ylim=c(0, 300)) + ggtitle("Muscle_Skeleta.p.funcExplorer") + xlab("-log10(p)") + theme(text=element_text(size=5))
fig.all = list(fig3, fig2, fig1, fig4)
ggsave("Muscle_Skeleta.p.funcExplorer.png", marrangeGrob(fig.all, nrow=1, ncol=length(fig.all), top = NULL), width = 10, height = 5)

p.obs[order(p.obs$V2), ][1:100, ]
p.obs2[order(p.obs2$V2), ][1:100, ]
p.obs3[order(p.obs3$V2), ][1:100, ]

p.obs[p.obs$V2<1e-8, ]
p.obs2[p.obs2$V2<1e-8, ]
p.obs3[p.obs3$V2<1e-8, ]

sort(setNames(sapply(sort(unique(p.obs[p.obs$V2<1e-8, V1])), function(x) paste0(strsplit(x, ":", fixed = T)[[1]][-1], collapse = ":")), NULL))
sort(setNames(sapply(sort(unique(p.obs2[p.obs2$V2<1e-8, V1])), function(x) paste0(strsplit(x, ":", fixed = T)[[1]][-1], collapse = ":")), NULL))


library(ggplot2)
require(reshape2)
library(ggpubr)

coexp.DGN = readRDS('coexp.DGN.rds')
coexp.Muscle = readRDS('coexp.Muscle.rds')

NmoduleMax.DGN = max(apply(coexp.DGN, 2, function(x) length(table(x))))
NmoduleMax.Muscle = max(apply(coexp.Muscle, 2, function(x) length(table(x))))

Dat.DGN = matrix(ncol = 6, nrow = NmoduleMax.DGN, dimnames = c(list(NULL, colnames(coexp.DGN) )))
Dat.DGN = apply(coexp.DGN, 2, function(x) {c(table(x), rep(NA, NmoduleMax.DGN-length(table(x))))})
DGN.0 = rbind(Dat.DGN[1, ], rep(colSums(Dat.DGN, na.rm = T)[1], 6), apply(coexp.DGN, 2, max)); rownames(DGN.0) = c("unclustered", "Total Gene", "NModule")
Dat.Muscle = matrix(ncol = 6, nrow = NmoduleMax.Muscle, dimnames = c(list(NULL, colnames(coexp.Muscle) )))
Dat.Muscle = apply(coexp.Muscle, 2, function(x) {c(table(x), rep(NA, NmoduleMax.Muscle-length(table(x))))})
Muscle.0 = rbind(Dat.Muscle[1, ], rep(colSums(Dat.Muscle, na.rm = T)[1], 6), apply(coexp.Muscle, 2, max)); rownames(Muscle.0) = c("unclustered", "Total Gene", "NModule")
Dat.all = rbind(cbind(as.data.frame(Dat.DGN[-1, ]), "Data" = rep("DGN", dim(Dat.DGN[-1, ])[1])), cbind(as.data.frame(Dat.Muscle[-1, ]), "Data" = rep("Muscle", dim(Dat.Muscle[-1, ])[1])) )

Dat.all = melt(Dat.all); colnames(Dat.all) = c("Data", "Method", "GeneinModule")
Dat.both = rbind(cbind(melt(DGN.0), "Data" = rep("DGN", dim(melt(DGN.0))[1])), cbind(melt(Muscle.0), "Data" = rep("Muscle", dim(melt(Muscle.0))[1])) )
colnames(Dat.both) = c("Var", "Method", "Num", "Data")

fig1 = ggplot(Dat.all, aes(x=Method, y=GeneinModule, col = Data)) + geom_boxplot() + coord_cartesian(ylim=c(0, 600)) + theme(axis.text=element_text(size=5))
fig3 = ggplot(Dat.both[Dat.both$Var != "NModule", ], aes(x = factor(Method), y = Num, col = Data, shape = Var, group = interaction(Data, Var))) + geom_line() + geom_point() + ggtitle("Num of Genes") + theme(axis.text=element_text(size=5))
fig4 = ggplot(Dat.both[Dat.both$Var == "NModule", ], aes(x = Method, y = Num, col = Data, group = Data)) + geom_line() + geom_point(shape = 7) + ggtitle("Num of Modules") + theme(axis.text=element_text(size=5)) + coord_cartesian(ylim=c(0, max(NmoduleMax.DGN, NmoduleMax.Muscle)))
ggsave("modules.png", ggarrange(fig1,                                                 # First row with scatter plot
                                ggarrange(fig3, fig4, ncol = 2, common.legend=T, labels = c("B", "C")), # Second row with box and dot plots
                                nrow = 2,
                                labels = "A")
       , width = 10, height = 5)






library(ggplot2)
library(gridExtra)
library(data.table)
require(dplyr)
require(reshape2)
library(ggpubr)

res.module.snp = matrix(nrow = 2, ncol = 6, dimnames = list(c("orig module", "smaller module"), c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5)))
res.uniq.snp = matrix(nrow = 2, ncol = 6, dimnames = list(c("orig module", "smaller module"), c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5)))
res.indpt.snp = matrix(nrow = 2, ncol = 6, dimnames = list(c("orig module", "smaller module"), c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5)))

file.q = paste0("/project2/xuanyao/llw/DGN_new/FDR/q.chr.module.perm", 1:10, '.rds')
p.obs = lapply(file.q, function(x) readRDS(x) )
for(fdr.level in c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5)){
  res = bind_rows(p.obs) %>% group_by(snp, p) %>% summarise("q" = mean(q)) %>% filter(q < fdr.level)
  module.snp = as.character(res$snp)
  snp = sort(unique(sapply(strsplit(module.snp, ":"), function(x) paste(x[-1], collapse = ":"))))

  res.module.snp["orig module", as.character(fdr.level)] = dim(res)[1]
  res.uniq.snp["orig module", as.character(fdr.level)] = length(snp)

  #write.table(snp, file = paste0("origmodule", fdr.level, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  res.indpt.snp["orig module", as.character(fdr.level)] = dim(read.table(paste0("indep.origmodule", fdr.level, ".txt")))[1]

}


file.q = paste0("FDR/q.chr.module.perm", 1:10, '.rds')
p.obs = lapply(file.q, function(x) readRDS(x) )
for(fdr.level in c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5)){
  res = bind_rows(p.obs) %>% group_by(snp, p) %>% summarise("q" = mean(q)) %>% filter(q < fdr.level)
  module.snp = as.character(res$snp)
  snp = sort(unique(sapply(strsplit(module.snp, ":"), function(x) paste(x[-1], collapse = ":"))))

  res.module.snp["smaller module", as.character(fdr.level)] = dim(res)[1]
  res.uniq.snp["smaller module", as.character(fdr.level)] = length(snp)

  #write.table(snp, file = paste0("smallermodule", fdr.level, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  res.indpt.snp["smaller module", as.character(fdr.level)] = dim(read.table(paste0("indep.smallermodule", fdr.level, ".txt")))[1]
}


res.module.snp
res.uniq.snp
res.indpt.snp

Dat = melt(res.module.snp); colnames(Dat) = c("module", "FDR", "Numofsignals")
fig1 = ggplot(Dat, aes(x=FDR, y=Numofsignals, col=module)) + geom_line() + geom_point() + coord_cartesian(ylim=c(0, 2000)) + ggtitle("(module snp) pair")

Dat = melt(res.uniq.snp); colnames(Dat) = c("module", "FDR", "Numofsignals")
fig2 = ggplot(Dat, aes(x=FDR, y=Numofsignals, col=module)) + geom_line() + geom_point() + coord_cartesian(ylim=c(0, 2000)) + ggtitle("unique snp")

Dat = melt(res.indpt.snp); colnames(Dat) = c("module", "FDR", "Numofsignals")
fig3 = ggplot(Dat, aes(x=FDR, y=Numofsignals, col=module)) + geom_line() + geom_point() + coord_cartesian(ylim=c(0, 300)) + ggtitle("independent snp")

fig.all = list(fig1, fig2, fig3)
ggsave("signalsv.s.FDR.png", ggarrange(plotlist=fig.all, nrow=1, ncol=3, common.legend=T), width = 10, height = 5)



{
  file1=$(ls origmodule0.*)

  for sig_uniq in $file1
  do
  echo $sig_uniq

  for chr in {1..22}
  do
  plink --bfile '/project2/xuanyao/llw/DGN/data/chr'$chr'_QCed' \
  --extract $sig_uniq \
  --indep-pairwise 50 5 0.2 \
  --out indep.eigene.eqtls.$chr
  echo chr$chr
  done

  cat indep.eigene.eqtls.*.prune.in > in.txt
  cat indep.eigene.eqtls.*.prune.out > out.txt

  grep -vFf out.txt $sig_uniq > 'indep.'$sig_uniq
  rm indep.eigene.eqtls* in.txt out.txt

  done

}


require(data.table)

meta = as.data.frame(fread('Muscle_Skeletal.gene_TSS.txt')); colnames(meta)[1]='chr'; rownames(meta) = meta$gene_id
signal = as.data.frame(fread('FDR/signals.chr.module.perm10.txt'))
cross_map = fread('/project2/xuanyao/data/GTEx_v8/mappability/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz')

signal$snp = sapply(signal$V1, function(x) paste(strsplit(x, ":", fixed = T)[[1]][-1], collapse = ":"))
signal$snpChr = sapply(signal$V1, function(x) paste0("chr", strsplit(x, ":", fixed = T)[[1]][2]))
signal$snpPos = as.numeric(sapply(signal$V1, function(x) strsplit(x, ":", fixed = T)[[1]][3]))
signal$module = sapply(signal$V1, function(x) strsplit(x, ":", fixed = T)[[1]][1])
signal$cisgene = apply(signal, 1, function(x){
  meta.tmp = meta[meta$chr == unlist(x['snpChr']), ]
  paste(meta.tmp[abs(meta.tmp[, 'start'] - as.numeric(x['snpPos']))<5*10^(5), 'gene_id'], collapse = ";")

})
signal$cisdis = apply(signal, 1, function(x){
  meta.tmp = meta[meta$chr == unlist(x['snpChr']), ]
  paste(abs(meta.tmp[, 'start'] - as.numeric(x['snpPos']))[abs(meta.tmp[, 'start'] - as.numeric(x['snpPos']))<5*10^(5)], collapse = ";")

})

signal$nearestgene = apply(signal, 1, function(x){
  meta.tmp = meta[meta$chr == x['snpChr'], ]
  meta.tmp[which.min(abs(meta.tmp[, 'start'] - as.numeric(x['snpPos']))), 'gene_id']

})
signal$nearestdis = apply(signal, 1, function(x){
  meta.tmp = meta[meta$chr == x['snpChr'], ]
  min(abs(meta.tmp[, 'start'] - as.numeric(x['snpPos'])))

})

coexp = readRDS('result/coexp.module.rds')
geneinModule = names(coexp$moduleLabels[coexp$moduleLabels==8])

num.crossinmodule = NULL
for(k in 1:20){
  a = strsplit(signal$cisgene[5], ";")[[1]][k]
  cross.gene = unique(cross_map[cross_map$V1 == a | cross_map$V2 == a, c(V1,V2)])
  num.crossinmodule = c(num.crossinmodule, sum(cross.gene %in% geneinModule[meta[geneinModule, 'chr']!='chr5']))
  #print(meta[meta$chr!='chr5'][cross.gene[cross.gene %in% geneinModule], ])
}

as.numeric(unlist(strsplit(signal$cisdis[5], ";")))[order(as.numeric(unlist(strsplit(signal$cisdis[5], ";"))))]
num.crossinmodule[order(as.numeric(unlist(strsplit(signal$cisdis[5], ";"))))]
meta[cross.gene[cross.gene %in% geneinModule[meta[geneinModule, 'chr']!='chr5']], ]

a = signal$nearestgene[5]
cross.gene = unique(cross_map[cross_map$V1 == a | cross_map$V2 == a, c(V1,V2)])
print(sum(cross.gene %in% geneinModule[meta[geneinModule, 'chr']!='chr5']))


cross_map.genes = unique(c(cross_map$V1, cross_map$V2))

res = data.frame(table(coexp$moduleLabels), stringsAsFactors = F); res=res[-1, ]; names(res) = c('module', "nGene")
res$nCross = rep(NA, 39)
for(k in 1:39){
  geneinModule = names(coexp$moduleLabels[coexp$moduleLabels==k])
  res[k, "nCross"] = sum(geneinModule %in% cross_map.genes)
}

res = rbind(res[, 1:2], res[, c(1,3)])
res$color = rep(c('All', 'nCross'), each = 39)
colnames(res) = c('module', 'nGene', 'color')
fig = ggplot(res,aes(module, nGene, fill=color))+
  geom_bar(stat="identity",position='dodge')
ggsave('test.png', fig)

res2
fig2 = ggplot(res2,aes(a, d))+
  geom_point() +
  geom_hline(yintercept=0.5, color='red') +
  ylim(0, 1) +
  labs(x = "module", y = "ratio")
ggsave('test.png', fig2)
  geom_abline(slope = 1, intercept = 0, color="red")

