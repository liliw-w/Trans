file.gene.annotation = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
file.mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_gene_mappability.txt.gz"
file.cross.mappability = "/project2/xuanyao/data/mappability/hg19_gencode19_75merExon_36merUTR_2mismatch_cross_mappability_symmetric_mean.txt.gz"

file.ex = "/project2/xuanyao/data/eQTLGen/DGN.rds"
file.covariates = "/project2/xuanyao/llw/DGN/data/covariates.txt"
file.ex.var.regressed = "ex.var.regressed.rds"

require(data.table)

## read mappability files
gene.annotation = fread(file.gene.annotation, header = TRUE)
low.mapp = fread(file.mappability)
cross.map = fread(file.cross.mappability)

gene.annotation$Geneid = sapply(gene.annotation$Geneid, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
low.mapp$V1 = sapply(low.mapp$V1, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
cross.map$V1 = sapply(cross.map$V1, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
cross.map$V2 = sapply(cross.map$V2, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])

## filter eQTLGen genes
gene.eQTLGen = fread("/project2/xuanyao/data/eQTLGen/eQTLGen.gene.txt", header = FALSE)

gnames = gene.eQTLGen$V1
gnames.meta = gnames[gnames %in% gene.annotation$Geneid]
rem.pseu = !(gene.annotation[match(gnames.meta, gene.annotation$Geneid), Class] %in% c("protein_coding", "lincRNA"))
rem.low.mapp = gnames.meta %in% low.mapp[low.mapp$V2 < 0.9, V1]
rem.cross.map = gnames.meta %in% unique(c(cross.map$V1, cross.map$V2))
rem.auto = !(gene.annotation[match(gnames.meta, gene.annotation$Geneid), Chromosome] %in% paste0("chr", 1:22))

GeneNameConv = gene.annotation[match(gnames.meta, gene.annotation$Geneid), GeneSymbol]
GeneName = gnames.meta
gnames.ensb = GeneName

res.eQTLGen = data.frame("gene" = GeneName,
                 "chr" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), Chromosome],
                 "start" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), Start],
                 "end" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), End],
                 "GeneNameConv" = GeneNameConv,
                 "Remove" = (rem.pseu | rem.low.mapp | rem.cross.map | rem.auto),
                 "RemovePseu" = rem.pseu,
                 "RemoveLowmapp" = rem.low.mapp,
                 "RemoveCrossmapp" = rem.cross.map,
                 "RemoveAuto" = rem.auto,
                 stringsAsFactors = FALSE)

write.table(as.matrix(res.eQTLGen),
            "eQTLGen.txt",
            sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

## filter DGN genes
ex = readRDS(file.ex)

gnames = colnames(ex)
gnames.meta = gnames[gnames %in% gene.annotation$GeneSymbol]
rem.pseu = !(gene.annotation[match(gnames.meta, gene.annotation$GeneSymbol), Class] %in% c("protein_coding", "lincRNA"))
rem.low.mapp = gnames.meta %in% gene.annotation[match(low.mapp[low.mapp$V2 < 0.9, V1], gene.annotation$Geneid), GeneSymbol]
rem.cross.map = gnames.meta %in% gene.annotation[match(unique(c(cross.map$V1, cross.map$V2)), gene.annotation$Geneid), GeneSymbol]
rem.auto = !(gene.annotation[match(gnames.meta, gene.annotation$GeneSymbol), Chromosome] %in% paste0("chr", 1:22))

GeneName = gnames[match(gnames.meta, gnames)]
GeneNameConv = gene.annotation[match(gnames.meta, gene.annotation$GeneSymbol), Geneid]
gnames.ensb = GeneNameConv

res.DGN = data.frame("gene" = GeneName,
                 "chr" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), Chromosome],
                 "start" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), Start],
                 "end" = gene.annotation[match(gnames.ensb, gene.annotation$Geneid), End],
                 "GeneNameConv" = GeneNameConv,
                 "Remove" = (rem.pseu | rem.low.mapp | rem.cross.map | rem.auto),
                 "RemovePseu" = rem.pseu,
                 "RemoveLowmapp" = rem.low.mapp,
                 "RemoveCrossmapp" = rem.cross.map,
                 "RemoveAuto" = rem.auto,
                 stringsAsFactors = FALSE)

write.table(as.matrix(res.DGN),
            "DGN.txt",
            sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


## Genes for final use
gnames.use = res.DGN$GeneNameConv[!res.DGN$Remove][res.DGN$GeneNameConv[!res.DGN$Remove] %in% res.eQTLGen$gene[!res.eQTLGen$Remove]]
gnames.use = res.DGN[match(gnames.use, res.DGN$GeneNameConv), "gene"]

## Part 2: regress out covariates
datExpr = ex[, gnames.use]

### read covariates
cov_all = t(as.matrix(read.table(file.covariates,
                                 sep = "\t", header = TRUE, row.names = 1,
                                 stringsAsFactors = FALSE, check.names = FALSE)))
### regress out covariates
extract_residual <- function(y, x){
  return(lm(y ~ x)$residuals)
}
ex_cov_regressed = apply(datExpr, 2, function(y) extract_residual(y, cov_all))

### save expression matrix with covariates regressed out
saveRDS(ex_cov_regressed, file.ex.var.regressed)


## Part 3: coexpressed modules
file.ex.var.regressed = "ex.var.regressed.rds"
file.coexp.module = 'coexp.module.rds'
minModuleSize = 20

library(WGCNA)


# load ex_cov_regressed
datExpr = readRDS(file.ex.var.regressed)

# Run WGCNA
### Parameter specification ###
minModuleSize = minModuleSize
MEDissThres = 0.15
if_plot_adjacency_mat_parameter_selection = F
if_plot_only_tree = F
if_plot_color_and_tree = F
if_plot_eigengene_heatmap_tree = F
if_heatmap_of_network = T


### Step1: network construction ###
# determine the paramter in adjacency function: pickSoftThreshold() OR pickHardThreshold()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower = sft$powerEstimate

# network construction
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


### Step2: module detection ###
# tree construction using hierarchical clustering based on TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

# branch cutting using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)

# eigene genes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

moduleLabels = match(mergedColors, c("grey", standardColors(100)))-1
names(moduleLabels) = colnames(datExpr)

print(table(moduleLabels))
cat("Number of modules:", max(moduleLabels), "\n")

# Save results
result = list(moduleColors = mergedColors,
              moduleLabels = moduleLabels,
              MEs = mergedMEs,
              old_moduleColors = dynamicColors,
              old_moduleLabels = dynamicMods,
              old_MEs = MEs,
              geneTree = geneTree)
saveRDS(result, file = file.coexp.module)



## split z matrices

library(data.table)

res.DGN = fread("/project2/xuanyao/data/eQTLGen/DGN.txt", header = TRUE)
res.DGN = as.data.frame(res.DGN)
result = readRDS("/scratch/midway2/liliw1/eQTLGen/result/coexp.module.rds")
Nmodule = max(result$moduleLabels)

for(chr in 1:22){
  z.chr.mat=fread(paste0("/scratch/midway2/liliw1/eQTLGen/z.chr", chr, ".txt.gz"), header=TRUE)
  z.chr.mat=as.matrix(z.chr.mat, rownames='snp')
  tmp = apply(z.chr.mat, 1, function(x) sum(is.na(x)))
  z.chr.mat = z.chr.mat[tmp == 0, ]
  
  for(m in 1:Nmodule){
    gene.in.module = res.DGN[match(names(result$moduleLabels)[result$moduleLabels == m], res.DGN$gene), "GeneNameConv"]
    z.module.chr.mat = cbind(data.frame("snp"=rownames(z.chr.mat)), z.chr.mat[, gene.in.module[gene.in.module %in% colnames(z.chr.mat)]])
    fwrite(z.module.chr.mat, paste0("/scratch/midway2/liliw1/eQTLGen/z/z.module", m, ".chr", chr, ".txt.gz"),
           sep = "\t", row.names = FALSE, col.names = TRUE)
    print(dim(z.module.chr.mat))
    cat("module:", m, "\n")
    
  }
  
  cat(chr, "\n")
}


## Part 3: pvalues
rm(list = ls())

file.ex.var.regressed = "/scratch/midway2/liliw1/eQTLGen/result/ex.var.regressed.rds"
file.coexp.module = "/scratch/midway2/liliw1/eQTLGen/result/coexp.module.rds"
params1 = "/scratch/midway2/liliw1/GTEx_v8/Whole_Blood/script/"
file.gene.meta = "/scratch/midway2/liliw1/eQTLGen/result/DGN.txt"

library(data.table)
datExpr = readRDS(file.ex.var.regressed)
gene.meta = read.table(file.gene.meta, sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
coexp.module = readRDS(file.coexp.module)
Nmodule = max(coexp.module$moduleLabels)

# Apply new PCO
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))


for(chr in 1:22){
  for(module in 1:Nmodule){
    file.z = paste0("/scratch/midway2/liliw1/eQTLGen/z/z.module", module, ".chr", chr, ".txt.gz")
    file.p = paste0("/scratch/midway2/liliw1/eQTLGen/p/p.module", module, ".chr", chr, ".rds")
    
    # zscore input
    z.mat = fread(file.z)
    z.mat = as.matrix(z.mat, rownames = TRUE)
    colnames(z.mat) = gene.meta[match(colnames(z.mat), gene.meta$GeneNameConv), "gene"]
    
    gene_in_cluster = data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
    gene_w_pos = merge(gene_in_cluster, gene.meta)
    gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]
    
    
    Sigma = cor(datExpr[, gene_trans])
    cat("Phase1 done.")
    SigmaO = ModifiedSigmaOEstimate(Sigma) 
    
    z.mat_trans = z.mat[, gene_trans]; rm(z.mat)
    p.all = ModifiedPCOMerged(Z.mat=z.mat_trans,
                              Sigma=Sigma, SigmaO=SigmaO)
    
    
    saveRDS(p.all, file = file.p)
    cat("module:", module, "\n")
  }
  
  cat(chr, "\n")
}


## Part 4: FDR

require(data.table)
require(qvalue)

file.signal = "/scratch/midway2/liliw1/eQTLGen/FDR/signals.qvalue.txt"

file.p = as.character(outer(1:19, 1:22, FUN = function(x, y) paste0("/scratch/midway2/liliw1/eQTLGen/p/p.module", x, ".chr", y, ".rds")))
p.obs = rbindlist(lapply(file.p, function(x) 
{tmp_y=readRDS(x);
as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)}))

pvalue_pc = p.obs$V2; names(pvalue_pc) = p.obs$V1
qobj_pc = qvalue(p = pvalue_pc, fdr.level=0.05)
print("FDR analysis of pc has been done!")

res = data.frame("snp" = names(qobj_pc$qvalues)[qobj_pc$significant],
                 "p" = qobj_pc$pvalues[qobj_pc$significant],
                 "q" = qobj_pc$qvalues[qobj_pc$significant],
                 stringsAsFactors = FALSE)

write.table(as.matrix(res),
            file.signal,
            sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



## Part 5: post analysis
rm(list = ls())

require(data.table)

file.signal = "/scratch/midway2/liliw1/eQTLGen/FDR/signals.qvalue.txt"
file.uniq = "/scratch/midway2/liliw1/eQTLGen/postanalysis/LD.prun.in.txt"
file.meta = "/project2/xuanyao/data/eQTLGen/meta.snp.txt.gz"

meta.snp = fread(file.meta, header = TRUE)
signals = read.table(file.signal,
                     sep = "\t", header = TRUE, row.names = 1,
                     stringsAsFactors = FALSE)
module.snp = rownames(signals)


snp = sort(unique(sapply(strsplit(module.snp, ":"), function(x) paste(x[-1], collapse = ":"))))
snp = sort(apply(meta.snp[match(snp, meta.snp$SNP), 2:3], 1, function(x) paste0(x, collapse = ":")))

write.table(snp, file = file.uniq, quote = FALSE, row.names = FALSE, col.names = FALSE)

dir_geno=/project2/xuanyao/data/GTEx_v8/genotype/
geno_prefix=chr
geno_suffix=_QCed
sig_uniq=postanalysis/LD.prun.in.txt
sig_indp=postanalysis/indep.signals.txt

module load plink
module load R/3.6.1

for chr in {1..22}
do
plink --bfile $dir_geno$geno_prefix$chr$geno_suffix \
--extract $sig_uniq \
--indep-pairwise 50 5 0.2 \
--out indep.eigene.eqtls.$chr
echo chr$chr
done

cat indep.eigene.eqtls.*.prune.in > in.txt
cat indep.eigene.eqtls.*.prune.out > out.txt

grep -vFf out.txt $sig_uniq > $sig_indp

rm indep.eigene.eqtls* in.txt out.txt

echo $dir_geno
echo "(module, snp):"$(wc $signals_file -l)
echo "unique snps:"$(wc $sig_uniq -l)
echo "independent snps:"$(wc $sig_indp -l)

