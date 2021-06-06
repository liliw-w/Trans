rm(list = ls())
parameter = commandArgs(trailingOnly = T)
PCO.script = parameter[1] # './script_lambda0.1/'
file.Sigma = parameter[2] # 'Sigma.DGN.module13_chr3.100.rds'
file.res = parameter[3] # "simulation.null.lambda0.1.K100.rds"

require(mvtnorm)
source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))

Sigma = as.matrix(readRDS(file.Sigma))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors

N.sample = 10^7
p.null.all  = list()

# z null
Z.null = rmvnorm(N.sample, rep(0, K), Sigma)

# PCO
p.null.all$'p.null.PCO' = as.numeric(ModifiedPCOMerged(Z.mat=Z.null, Sigma=Sigma, SigmaO=SigmaO))

cat("PCO done.")
saveRDS(p.null.all, file = file.res)

# PC1
PC1 = Z.null %*% eigen.vec[, 1]
p.null.all$'p.null.PC1' = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

cat("PC1 done.")

# univariate minp
p.null.all$'p.null.minp' = apply(Z.null, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

cat("minp done.")
saveRDS(p.null.all, file = file.res)
