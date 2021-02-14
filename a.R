library(data.table)

chr=9
dat=fread(paste0("chr",chr,".geno.raw"), header=T, sep = " ", quote = FALSE)
geno=as.matrix(dat[,7:ncol(dat)])
rownames(geno)=dat[, IID]
colnames(geno)=unlist(lapply(strsplit(colnames(dat)[7:ncol(dat)],"_"),"[",1))
fwrite(geno,paste0("chr",chr,".genotype.matrix.eqtl.txt"), sep="\t", row.names = T, col.names = T)

G = fread(paste0("chr",chr,".genotype.matrix.eqtl.txt"))

Y = readRDS("result/ex.var.regressed.rds")
samples = rownames(Y)

snp = "9:102349016"; gene = "ENSG00000002745.12"
y = Y[samples, gene]
x = G[samples, snp]
summary(lm(y~x))$coefficient[2, 3]


params1 = "script/"
source(paste0(params1, "ModifiedPCOMerged.R"))
source(paste0(params1, "liu.R"))
source(paste0(params1, "liumod.R"))
source(paste0(params1, "davies.R"))
dyn.load(paste0(params1, "qfc.so"))
source(paste0(params1, "ModifiedSigmaOEstimate.R"))


library(data.table)

module = 12; chr = 3
p = readRDS("p/p.module12.chr3.rds")
zmat = fread("z/z.module12.chr3.txt.gz", header = TRUE)

zmat = as.matrix(zmat, rownames = 'snp')

datExpr = readRDS("result/ex.var.regressed.rds")
coexp.module = readRDS("result/coexp.module.rds")
gene.meta = read.table("result/gene.meta.txt", sep = '\t', header = T, row.names = NULL, stringsAsFactors = F,  check.names = FALSE)
gene_in_cluster = data.frame("gene" = names(coexp.module$moduleLabels[coexp.module$moduleLabels == module]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene.meta)
gene_trans = gene_w_pos[gene_w_pos[, "chr"] != paste0("chr", chr), "gene"]


Sigma = cor(datExpr[, gene_trans])
SigmaO = ModifiedSigmaOEstimate(Sigma)

tmp = c("3:56815721", "3:100001720")
z.mat_trans = zmat[tmp, gene_trans]
p.all = ModifiedPCOMerged(Z.mat=z.mat_trans,
                          Sigma=Sigma, SigmaO=SigmaO)
p.all

z.mat_trans = zmat[1:100, gene_trans]
p.all = ModifiedPCOMerged(Z.mat=z.mat_trans,
                          Sigma=Sigma, SigmaO=SigmaO)

cbind(p.all, p[names(p.all)])[1:50, ]

q1 = readRDS("FDR/q.chr.module.perm1.rds")
q1[order(q1$p), ][1:10, ]

tmp = as.character(q1[order(q1$p), 'snp'])[2:4]


library(ggplot2)
library(gridExtra)
library(data.table)


Nmodule = 32; plot.name = "plots/Muscle_Skeletal_WGCNA.p.png"

file.p = as.character(outer(1:Nmodule, 1:22, FUN = function(x, y) paste0("p/p.module", x, ".chr", y, ".rds")))
p.obs = rbindlist(lapply(file.p, function(x)
{tmp_y=readRDS(x);
tmp_y = tmp_y[tmp_y<1e-6];
if(length(tmp_y > 0)){as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)}
}))

file.q = as.character(outer(1:Nmodule, 1:22, FUN = function(x, y) paste0("FDR/q.module", x, ".chr", y, ".perm3.rds")))
p.obs = unlist(rbind(lapply(file.q, function(x)
{tmp_y=readRDS(x); setNames(tmp_y$p, tmp_y$snp)})))
p.obs = data.frame("V1" = names(p.obs), "V2" = p.obs, stringsAsFactors = FALSE)

fsort(p.obs$V2)[1:50]

p.obs$V2[p.obs$V2==0] = 1e-20
fig1 = ggplot(p.obs[p.obs$V2<1e-6, ], aes(x=-log10(V2))) + geom_histogram(binwidth=0.5) + coord_cartesian(xlim=c(0,20), ylim=c(0, 250))
fig2 = ggplot(p.obs[p.obs$V2<1e-6, ], aes(y=-log10(V2))) + geom_boxplot()+ coord_cartesian(ylim=c(0,20))
fig.all = list(fig1, fig2)
ggsave(plot.name, marrangeGrob(fig.all, nrow=1, ncol=2, top = NULL), width = 10, height = 5)
sum(p.obs$V2<1e-6)
