# Z alternative
# FDR based on empirical pnull

rm(list = ls())
parameter = commandArgs(trailingOnly = T)

change = parameter[1] #'N'
oracle.thre = as.numeric(parameter[2]) #0.1
PCO.script = parameter[3] #'./script_lambda0.1/'
file.Sigma = parameter[4] #'Sigma.DGN.module13_chr3.100.rds'
file.res = parameter[5] #'simulation.alt.N.lambda0.1.K100.rds'
var.b = as.numeric(parameter[6]) #0.003
caus = as.numeric(parameter[7]) #0.3
N = as.numeric(parameter[8]) #500

require(mvtnorm)

source(paste0(PCO.script, "ModifiedPCOMerged.R"))
source(paste0(PCO.script, "liu.R"))
source(paste0(PCO.script, "liumod.R"))
source(paste0(PCO.script, "davies.R"))
dyn.load(paste0(PCO.script, "qfc.so"))
source(paste0(PCO.script, "ModifiedSigmaOEstimate.R"))


N.sample = 10^4; N.sim = 10^3
#thre = 0.05/(50*1000000)

Sigma = as.matrix(readRDS(file.Sigma))
SigmaO = ModifiedSigmaOEstimate(Sigma)
K = dim(Sigma)[1]
eigen.res = eigen(Sigma)
lambdas = eigen.res$values
eigen.vec = eigen.res$vectors
Sigma.trunc.inv = eigen.vec[, which(lambdas>oracle.thre)] %*% diag(1/lambdas[lambdas>oracle.thre]) %*% t(eigen.vec[, which(lambdas>oracle.thre)])

if(change == "N"){
  cat(paste0("Various ", change), '\n')

  # Different true effects
  N.seq = c(200, 400, 600, 800)
  models = paste0("N=", N.seq)
  #c("u1", "uPCO", "uK", paste0("N=", N.seq), "100%", "70%", "30%")


  res.alt = NULL
  p_alt_all = array(dim = c(N.sample, 4, length(models), N.sim),
                    dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    #beta.n = as.matrix(rnorm(K, sd = sqrt(var.b)))
    beta.n = rbind(as.matrix(rnorm(floor(caus*K), sd = sqrt(var.b))), as.matrix(rep(0, K-floor(caus*K))))
    B[, models] = beta.n %*% sqrt(N.seq)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% Sigma.trunc.inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] = 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% Sigma.trunc.inv %*% B[, model])))^2, 1)

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p_alt_all[, method.tmp, model, i] = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

      # univariate minp
      method.tmp = "minp"
      p_alt_all[, method.tmp, model, i] = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

      # PCO
      method.tmp = "PCO"
      p_alt_all[, method.tmp, model, i] = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file.res)}
  }

}else if(change == "caus"){
  cat(paste0("Various ", change), '\n')

  caus.seq = c(1, 5, 10, 30, 50)
  models = paste0("caus=", caus.seq)

  res.alt = NULL
  p_alt_all = array(dim = c(N.sample, 4, length(models), N.sim),
                    dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    B[, models] = sapply(caus.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b)), rep(0, K-x))) ) * sqrt(N)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% Sigma.trunc.inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] = 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% Sigma.trunc.inv %*% B[, model])))^2, 1)

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p_alt_all[, method.tmp, model, i] = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

      # univariate minp
      method.tmp = "minp"
      p_alt_all[, method.tmp, model, i] = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

      # PCO
      method.tmp = "PCO"
      p_alt_all[, method.tmp, model, i] = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))

    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file.res)}
  }

}else if(change == "caus_fix"){
  cat(paste0("Various ", change), '\n')

  caus.seq = c(1, 5, 10, 30, 50)
  models = paste0("caus.fix=", caus.seq)

  res.alt = NULL
  p_alt_all = array(dim = c(N.sample, 4, length(models), N.sim),
                    dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    B[, models] = sapply(caus.seq, function(x) as.matrix(c(rnorm(x, sd = sqrt(var.b/x)), rep(0, K-x))) ) * sqrt(N)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% Sigma.trunc.inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] = 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% Sigma.trunc.inv %*% B[, model])))^2, 1)

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p_alt_all[, method.tmp, model, i] = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

      # univariate minp
      method.tmp = "minp"
      p_alt_all[, method.tmp, model, i] = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

      # PCO
      method.tmp = "PCO"
      p_alt_all[, method.tmp, model, i] = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))

    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file.res)}
  }

}else if(change == "var"){
  cat(paste0("Various ", change), '\n')

  # Different true effects
  var.seq = c(0.005, 0.01, 0.05, 0.1, 0.2)
  models = paste0("var=", var.seq)

  res.alt = NULL
  p_alt_all = array(dim = c(N.sample, 4, length(models), N.sim),
                    dimnames = list(NULL, c("Oracle", "PC1", "minp", "PCO"), models, NULL))
  for(i in 1:N.sim){
    B = matrix(rep(NA, K*length(models)), ncol = length(models),
               dimnames = list(NULL, models))
    B[, models] = sapply(var.seq, function(x)
      as.matrix(c(rnorm(floor(caus*K), sd = sqrt(x)), rep(0, K-floor(caus*K)) ) ) ) * sqrt(N)

    for(model in models){
      Z.alt = rmvnorm(N.sample, B[, model], Sigma)

      # Oracle
      method.tmp = "Oracle"
      T.oracle = as.numeric(Z.alt %*% Sigma.trunc.inv %*% B[, model])
      p_alt_all[, method.tmp, model, i] = 1-pchisq((T.oracle/as.numeric(sqrt(t(B)[model, ] %*% Sigma.trunc.inv %*% B[, model])))^2, 1)

      # PC1
      method.tmp = "PC1"
      PC1 = Z.alt %*% eigen.vec[, 1]
      p_alt_all[, method.tmp, model, i] = as.numeric(2*pnorm(-abs(PC1/sqrt(lambdas[1]))))

      # univariate minp
      method.tmp = "minp"
      p_alt_all[, method.tmp, model, i] = apply(Z.alt, 1, function(x) min(1-pchisq(x^2, 1))*length(x) )

      # PCO
      method.tmp = "PCO"
      p_alt_all[, method.tmp, model, i] = as.numeric(ModifiedPCOMerged(Z.mat=Z.alt, Sigma=Sigma, SigmaO=SigmaO))
    }
    cat("Simulation: ", i, "\n")
    if(i %% 20 == 0){
      saveRDS(p_alt_all, file.res)}
  }

}else{cat("Wrong parameter!")}
