rm(list = ls())
PCO.script = './script_lambda0.1/'
file_Sigma = './new_Sigma/Sigma-new_DGN_module29_K101.rds'
file_dat_null = "./new_Sigma/simulation.null.lambda0.1.K101.rds"


library(mvtnorm)
source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))

Sigma = as.matrix(readRDS(file_Sigma))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors

N.sample = 10^7
p_null_all  = list()

# z null
Z.null = rmvnorm(N.sample, rep(0, K), Sigma)

# PCO
p_null_all$'p.null.PCO' = as.numeric(ModifiedPCOMerged(Z.mat=Z.null, Sigma=Sigma, SigmaO=SigmaO))

cat("PCO done.")
saveRDS(p_null_all, file = file_dat_null)

# PC1
PC1 = Z.null %*% eigen.vec[, 1]
p_null_all$'p.null.PC1' = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

cat("PC1 done.")
saveRDS(p_null_all, file = file_dat_null)

# univariate minp
p_null_all$'p.null.minp' = apply(Z.null, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

cat("minp done.")
saveRDS(p_null_all, file = file_dat_null)

